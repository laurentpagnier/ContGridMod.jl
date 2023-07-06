export find_gen, find_node, disc_dynamics

"""
    disc_dynamics(dm::DiscModel, tstart::Real, tend::Real, delta_p::Union{Real,Array{Real,1}}[, faultid::Int=0, coords::Array{<:Real,2}=[NaN NaN], dt::Real=1e-2, scale_factor::Real=0.0, tol::Real=1e-10, maxiter::Int=30, dmin::Real=1e-4, alg=TRBDF2(), solve_kwargs::Dict{String,Any}=Dict{String,Any}()])::ODESolution

This function allows to either specify the GPS coordinates of where the fault occurs or the ID of the generator.
If both are given the ID will be used. If the GPS coordinates are given, the scale_factor must be given.
Alternatively, an array of size dm.Nbus can be provided which specifies the change in power for each node.
"""
function disc_dynamics(
    dm::DiscModel,
    tstart::Real,
    tend::Real,
    delta_p::Union{Real,Array{Real,1}};
    faultid::Int=0,
    coords::Array{<:Real,2}=[NaN NaN],
    dt::Real=1e-2,  # Distance of time steps at which the solution is returned
    scale_factor::Real=0.0,
    tol::Real=1e-10,  # Target tolerance for the Newtoon-Raphson solver
    maxiter::Int=30,  # Maximum iteration for the Newton-Raphson solver
    dmin::Real=1e-4,  # Minimum amount of damping for load buses
    alg::OrdinaryDiffEqAlgorithm=TRBDF2(),  # The solver that is passed to the solve function
    solve_kwargs::Dict{String,Any}=Dict{String,Any}()  # Keyword arguments to the solve function
)::ODESolution
    if (faultid == 0 && size(delta_p, 1) != dm.Nbus)
        if (all(isnan.(coords)) || scale_factor == 0.0)
            throw(ErrorException("Either the coordinates or the ID of the faulty generator need to be specified. If the coordinates are given, scale_factor must be specified."))
        end
        _, tmp = find_gen(dm, coords, delta_p, scale_factor=scale_factor)
        faultid = tmp[1]
    end

    # Data preperation
    Nbus = dm.Nbus
    is_producing = dm.p_gen .> 1e-4
    id_gen = sort(dm.id_gen[is_producing])
    id_load = sort(setdiff(1:Nbus, id_gen))
    ng = length(id_gen)
    p_gen = zeros(Nbus)
    p_gen[dm.id_gen] .= dm.p_gen
    p = p_gen - dm.p_load
    p .-= sum(p) / Nbus
    nline = dm.Nline
    inc = sparse([dm.id_line[:, 1]; dm.id_line[:, 2]], [1:nline; 1:nline], [-ones(nline); ones(nline)])

    # create arrays of dynamical parameters
    mg = dm.m_gen[is_producing]
    dg = dm.d_gen[is_producing] + dm.d_load[id_gen]
    dl = max.(dm.d_load[id_load], dmin)

    # get the stable solution
    Ybus = -im * inc * sparse(1:nline, 1:nline, dm.b) * inc'
    q = zeros(Nbus)
    V = ones(Nbus)
    theta = zeros(Nbus)
    V, theta, _ = NRsolver(Ybus, V, theta, p, q, Array{Int64,1}([]), dm.id_slack, tol=tol, maxiter=maxiter)

    # preparations for dynamical simulation, mass matrix and fault vectors
    if (faultid != 0)
        dp = zeros(Nbus)
        dp[faultid] = delta_p
    else
        dp = delta_p
    end
    ix = collect(1:Nbus)
    for (i, j) in enumerate(id_gen)
        ix[j] = Nbus + i
    end
    mass = ones(Nbus)
    mass[id_load] = dl
    mass = 1 ./ [mass; mg]
    u0 = [theta; zeros(ng)]
    tspan = (tstart, tend)

    function swing!(du, u, _, _)
        du[id_gen] = u[Nbus+1:end]
        du[ix] = p + dp - inc * ((dm.b .* sin.(inc' * u[1:Nbus])))
        du[Nbus+1:end] -= dg .* u[Nbus+1:end]
        du .*= mass
    end

    # Easy indexing for Jacobian calculations
    id1 = [id_load; Nbus+1:Nbus+ng]
    id2 = [id_gen; id_load]
    id3 = [id_load; id_gen]
    function jacobian!(J, u, _, _)
        J0 = -inc * ((dm.b .* cos.(inc' * u[1:Nbus])) .* inc')
        J[:, :] = spzeros(Nbus + ng, Nbus + ng)
        J[id1, id2] = J0[id3, id2]
        J[Nbus+1:end, Nbus+1:end] = -spdiagm(dg)
        J[:, :] = mass .* J
    end

    jac_proto = spzeros(Nbus + ng, Nbus + ng)
    jacobian!(jac_proto, ones(Nbus + ng), 0, 0)
    for i = 1:Nbus+ng
        jac_proto[i, i] += 1
    end
    # save the frequencies at the predefined time steps
    tt = tstart:dt:tend

    # solve the swing equations
    func = ODEFunction(swing!, jac=jacobian!, jac_prototype=jac_proto)
    prob = ODEProblem(func, u0, tspan)
    sol = solve(prob, alg, tstops=tt, solve_kwargs...)
    return sol
end

"""
    
    NRsolver(Ybus::SparseMatrixCSC{<:Complex,<:Int}, V::Array{<:Real,1}, theta::Array{<:Real,1}, p::Array{<:Real,1}, q::Array{<:Real,1}, idpq::Array{<:Int,1}, id_slack::Int; tol::Real=1E-6, maxiter::Int=14)::Tuple{Array{<:Real,1},Array{<:Real,1},Int}

Use the Newton Raphson method to solve the powerflow equations.This method is
adapted from its version on the Pantagruel repository (https://doi.org/10.5281/zenodo.2642175).
For information on solving the power flow equations with Newton-Raphson, see, 
for instance, V. Vittal and A. Bergen, Power systems analysis, Prentice Hall, 1999. 
"""
function NRsolver(
    Ybus::SparseMatrixCSC{<:Complex,<:Int},
    V::Array{<:Real,1},
    theta::Array{<:Real,1},
    p::Array{<:Real,1},
    q::Array{<:Real,1},
    idpq::Array{<:Int,1},
    id_slack::Int;
    tol::Real=1E-6,
    maxiter::Int=14
)::Tuple{Array{<:Real,1},Array{<:Real,1},Int}
    #=
    This method is adapted from its version on the
    Pantagruel repository (https://doi.org/10.5281/zenodo.2642175).
    For information on solving the power flow equations with 
    Newton-Raphson, see, for instance,
    V. Vittal and A. Bergen, Power systems analysis,
    Prentice Hall, 1999.
    =#
    nb = size(Ybus, 1)
    error = 2 * tol
    iter = 0
    id = [1:id_slack-1; id_slack+1:nb]
    while (error > tol && iter < maxiter)
        Vc = V .* exp.(im * theta)
        S = Vc .* conj(Ybus * Vc)
        dPQ = [real(S[id]) - p[id]; imag(S[idpq]) - q[idpq]]
        dsdth = -im * sparse(1:nb, 1:nb, Vc) * conj(Ybus) * sparse(1:nb, 1:nb, conj(Vc)) +
                im * sparse(1:nb, 1:nb, Vc .* conj(Ybus * Vc))
        dsdv = sparse(1:nb, 1:nb, Vc) * conj(Ybus) * sparse(1:nb, 1:nb, exp.(-im * theta)) +
               sparse(1:nb, 1:nb, exp.(im * theta) .* conj(Ybus * Vc))
        J = [real(dsdth[id, id]) real(dsdv[id, idpq]); imag(dsdth[idpq, id]) imag(dsdv[idpq, idpq])]
        x = J \ dPQ
        theta[id] = theta[id] - x[1:nb-1]
        if (!isempty(idpq))
            V[idpq] -= x[nb:end]
        end
        error = maximum(abs.(dPQ))
        iter += 1
    end
    if (iter == maxiter)
        println("Max iteration reached, error: ", error)
    end
    return V, theta, iter
end

"""
    find_node(grid_coord::Array{<:Real,2}, disc_coord::Array{<:Real,2})

Find the the nodes that are closest to given GPS coordinates.
"""
function find_node(
    grid_coord::Array{<:Real,2},
    disc_coord::Array{<:Real,2}
)::Array{<:Int,2}
    ids = Int64.(zeros(size(disc_coord, 1)))
    for i in eachindex(ids[:, 1])
        ids[i] = argmin((grid_coord[:, 1] .- disc_coord[i, 1]) .^ 2 +
                        (grid_coord[:, 2] .- disc_coord[i, 2]) .^ 2)
    end

    return ids
end

"""
    find_gen(dm::DiscModel, gps_coord::Array{<:Real,2}, dP::Real[, scale_factor::Real=1.0])::Array{<:Int,2}

Find the generators that are closest to given GPS coordinates and that produce more than
dP.
The GPS coordinate is first projected using the Alber projection and then scaled by the given
scale_factor.
"""
function find_gen(
    dm::DiscModel,
    gps_coord::Array{<:Real,2},
    dP::Real;
    scale_factor::Real=1.0
)::Array{<:Int,2}
    #find the the nearest generator that can "withstand" a dP fault
    coord = albers_projection(gps_coord[:, [2; 1]] ./ (180 / pi)) / scale_factor
    idprod = findall((dm.p_gen .> 0.0))
    idavail = findall(dm.max_gen[idprod] .> abs(dP)) # find id large gens in the prod list
    id = Int64.(zeros(size(coord, 1)))  # index in the full gen list
    for i in eachindex(coord[:, 1])
        temp = dm.id_gen[idprod[idavail]]
        id[i] = idprod[idavail[argmin((dm.coord[temp, 1] .- coord[i, 1]) .^ 2 +
                                      (dm.coord[temp, 2] .- coord[i, 2]) .^ 2)]]
    end
    return dm.id_gen[id, :]
end
