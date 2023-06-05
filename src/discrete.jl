using LinearAlgebra
using SparseArrays
using IterativeSolvers
using DifferentialEquations

export find_gen, find_node, find_node2, NRsolver, disc_dynamics

function disc_dynamics(
    dm::DiscModel,
    tstart::Float64,
    tend::Float64,
    delta_p::Union{Float64,Array{Float64,1}};
    faultid::Int64=0,
    coords::Array{Float64,2}=[NaN NaN],
    dt::Float64=1e-2,  # Distance of time steps at which the solution is returned
    scale_factor::Float64=0.0,
    tol::Float64=1e-10,  # Target tolerance for the Newtoon-Raphson solver
    maxiter::Int64=30,  # Maximum iteration for the Newton-Raphson solver
    dmin::Float64=1e-4,  # Minimum amount of damping for load buses
    reorder::Bool=false,  # Whether the results will be reorder such that the generators are first 
    alg=TRBDF2(),  # The solver that is passed to the solve function
    solve_kwargs::Dict=Dict()  # Keyword arguments to the solve function
)
    #=
    This function allows to either specify the GPS coordinates of where the fault occurs or the ID of the generator.
    If both are given the ID will be used. If the GPS coordinates are given, the scale_factor must be given.
    Alternatively, an array of size dm.Nbus can be provided which specifies the change in power for each node.
    =#
    if (faultid == 0 && size(delta_p, 1) != dm.Nbus)
        if (coords == [NaN NaN] || scale_factor == 0.0)
            throw(ErrorException("Either the coordinates or the ID of the faulty generator need to be specified. If the coordinates are given, scale_factor must be specified."))
        end
        ~, tmp = find_gen(dm, coords, dP, scale_factor=scale_factor)
        faultid = tmp[1]
    end

    # Data preperation
    Nbus = dm.Nbus
    is_producing = dm.p_gen .> 0
    id_gen = dm.id_gen[is_producing]
    id_load = setdiff(1:Nbus, id_gen)
    ng = length(id_gen)
    nl = length(id_load)
    p_gen = zeros(Nbus)
    p_gen[dm.id_gen] .= dm.p_gen
    p = p_gen - dm.p_load
    p .-= sum(p) / Nbus
    edges = Int64.(zeros(size(dm.id_line)))
    line_start = dm.id_line[:, 1]
    line_end = dm.id_line[:, 2]

    # bus reordering: generator buses first 
    for i in 1:ng
        edges[line_start.==id_gen[i], 1] .= i
        edges[line_end.==id_gen[i], 2] .= i
    end
    for i in 1:nl
        edges[line_start.==id_load[i], 1] .= i + ng
        edges[line_end.==id_load[i], 2] .= i + ng
    end
    nline = dm.Nline
    inc = sparse([edges[:, 1]; edges[:, 2]], [1:nline; 1:nline], [-ones(nline); ones(nline)])

    # create arrays of dynamical parameters
    mg = dm.m_gen[is_producing]
    dg = dm.d_gen[is_producing] + dm.d_load[id_gen]
    dl = max.(dm.d_load[id_load], dmin)
    pg = p[id_gen]
    pl = p[id_load]
    p = [pg; pl]

    # get the stable solution
    Ybus = -im * inc * sparse(1:nline, 1:nline, dm.b) * inc'
    q = zeros(Nbus)
    V = ones(Nbus)
    theta = zeros(Nbus)
    V, theta, ~ = NRsolver(Ybus, V, theta, p, q, Array{Int64,1}([]), 1, tol=tol, maxiter=maxiter)

    # preparations for dynamical simulation
    if (faultid != 0)
        dp = zeros(Nbus)
        dp[faultid] = delta_p
    else
        dp = delta_p
    end
    ix = [Nbus+1:Nbus+ng; ng+1:Nbus]
    mass = 1 ./ [ones(ng); dl; mg]
    u0 = [theta; zeros(ng)]
    tspan = (tstart, tend)

    function swing!(du, u, para, t)
        du[1:ng] = u[Nbus+1:end]
        du[ix] = p + dp - inc * ((dm.b .* sin.(inc' * u[1:Nbus])))
        du[Nbus+1:end] -= dg .* u[Nbus+1:end]
        du .*= mass
    end

    function jacobian!(J, u, p, t)
        J0 = -inc * ((dm.b .* cos.(inc' * u[1:Nbus])) .* inc')
        J[:, :] = [spzeros(ng, Nbus) sparse(1.0I, ng, ng); J0[ng+1:end, 1:ng] J0[ng+1:Nbus, ng+1:Nbus] spzeros(nl, ng); J0[1:ng, 1:ng] J0[1:ng, ng+1:Nbus] -dg.*sparse(1.0I, ng, ng)]
        J[:, :] = mass .* J
    end
    jac_proto = spzeros(Nbus + ng, Nbus + ng)
    jacobian!(jac_proto, ones(Nbus + ng), 0, 0)
    for i=1:Nbus + ng
        jac_proto[i, i] += 1
    end
    # save the frequencies at the predefined time steps
    saved_values = SavedValues(Float64, Vector{Float64})
    tt = tstart:dt:tend
    cb = SavingCallback((u, t, integrator) -> [u[Nbus+1:end]; integrator(t, Val{1}, idxs=ng+1:Nbus)], saved_values, saveat=tt)

    # solve the swing equations
    func = ODEFunction(swing!, jac=jacobian!, jac_prototype=jac_proto)
    prob = ODEProblem(func, u0, tspan)
    sol = solve(prob, alg, tstops=tt, callback=cb; solve_kwargs...)
    vals = reduce(hcat, saved_values.saveval)'

    # get ids to reorder for return values 
    re_id = Array(1:Nbus)
    if !reorder
        for (i, j) in enumerate(id_gen)
            re_id[j] = i
        end
        for (i, j) in enumerate(id_load)
            re_id[j] = i + ng
        end
    end
    return sol
    # return saved_values.t, vals[:, re_id]
end


function NRsolver(
    Ybus::SparseMatrixCSC{ComplexF64,Int64},
    V::Array{Float64,1},
    theta::Array{Float64,1},
    p::Array{Float64,1},
    q::Array{Float64,1},
    idpq::Array{Int64,1},
    id_slack::Int64;
    tol::Float64=1E-6,
    maxiter::Int64=14
)
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
            V[idpq] -= x[n:end]
        end
        error = maximum(abs.(dPQ))
        iter += 1
    end
    if (iter == maxiter)
        println("Max iteration reached, error: ", error)
    end
    return V, theta, iter
end



function get_discrete_id(
    grid_coord::Array{Float64,2},
    disc_coord::Array{Float64,2},
)
    ids = Int64.(zeros(size(disc_coord, 1)))
    for i in eachindex(ids[:, 1])
        ids[i] = argmin((grid_coord[:, 1] .- disc_coord[i, 1]) .^ 2 +
                        (grid_coord[:, 2] .- disc_coord[i, 2]) .^ 2)
    end

    return ids
end


function get_discrete_values(
    ids::Vector{Int64},
    v::Vector{Float64}
)
    return [v[ids[i]] for i = 1:length(ids)]
end


function get_discrete_values(
    grid_coord::Matrix{Float64},
    disc_coord::Matrix{Float64},
    v::Vector{Float64}
)
    disc_v = zeros(size(disc_coord, 1))
    for i in eachindex(disc_v[:, 1])
        id = argmin((grid_coord[:, 1] .- disc_coord[i, 1]) .^ 2 +
                    (grid_coord[:, 2] .- disc_coord[i, 2]) .^ 2)
        disc_v[i] = disc_v[i] + v[id] # changed to avoid "Mutating arrays is not supported"
    end
    return disc_v
end


function find_gen(
    dm::DiscModel,
    gps_coord::Array{Float64,2},
    dP::Float64;
    scale_factor::Float64=1.0
)
    #find the the nearest generator that can "withstand" a dP fault
    coord = albers_projection(gps_coord[:, [2; 1]] ./ (180 / pi))
    coord = coord[:, [2, 1]] / scale_factor
    idprod = findall((dm.p_gen .> 0.0))
    idavail = findall(dm.max_gen[idprod] .> dP) # find id large gens in the prod list
    # println(size(idprod))
    # println(idavail)
    #idprod = findall((dm.gen .> 0.0))
    id = Int64.(zeros(size(coord, 1)))  # index in the full gen list
    id2 = Int64.(zeros(size(coord, 1))) # index in the producing gen list
    for i in eachindex(coord[:, 1])
        temp = dm.id_gen[idprod[idavail]]
        id[i] = idprod[idavail[argmin((dm.coord[temp, 1] .- coord[i, 1]) .^ 2 +
                                      (dm.coord[temp, 2] .- coord[i, 2]) .^ 2)]]
        id2[i] = idavail[argmin((dm.coord[temp, 1] .- coord[i, 1]) .^ 2 +
                                (dm.coord[temp, 2] .- coord[i, 2]) .^ 2)]
    end
    return id, id2
end
