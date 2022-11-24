using LinearAlgebra
using SparseArrays
using IterativeSolvers
using DifferentialEquations

export find_gen, find_node, find_node2, NRsolver, disc_dynamics

function disc_dynamics(
    dm::DiscModel,
    tstart::Number,
    tend::Number,
    delta_p::Union{Float64, Vector{Float64}};
    faultid::Int64 = 0,
    coords::Matrix{Float64} = [NaN NaN],
    dt::Float64 = 1e-2,  # Distance of time steps at which the solution is returned
    scale_factor::Float64 = 0.0,
    tol::Float64 = 1e-10,  # Target tolerance for the Newtoon-Raphson solver
    maxiter::Int64 = 30,  # Maximum iteration for the Newton-Raphson solver
    alg=TRBDF2(),  # The solver that is passed to the solve function
    solve_kwargs::Dict = Dict()  # Keyword arguments to the solve function
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


    nline = dm.Nline
    inc = sparse([dm.id_line[:, 1]; dm.id_line[:, 2]], [1:nline; 1:nline], [-ones(nline); ones(nline)])

    # create arrays of dynamical parameters
    mg = dm.m_gen[is_producing]
    dg = dm.d_gen[is_producing] + dm.d_load[id_gen]
    dl = max.(dm.d_load[id_load], dmin)

    # get the stable solution
    theta, _ = compute_phases(dm.b, dm.id_line, pref, theta = dm.th,
        id_slack = 1, tol = 1E-6, maxiter = 14)

    # preparations for dynamical simulation
    if (faultid != 0)
        dp = zeros(Nbus)
        dp[faultid] = delta_p
    else
        dp = delta_p
    end

    temp = ones(Nbus)
    temp[id_load] = dl
    mass = [temp; mg]
    u0 = [theta; zeros(ng)]
    tspan = (tstart, tend)
    function swing!(du, u, param, t)
        # u[[1:Nbus] phases
        # u[Nbus+1:end] generator freqs.
        dpe = p + dp + inc * ((dm.b .* sin.(inc' * u[1:Nbus])))
        du[id_gen] = u[Nbus+1:end]
        du[id_load] = dpe[id_load]
        du[Nbus+1:end] = dpe[id_gen] - dg .* u[Nbus+1:end]
        du ./= mass
    end

    function jacobian!(J, u, param, t)
        J0 = inc * ((dm.b .* cos.(inc' * u[1:Nbus])) .* inc')
        J[id_load, :] = [J0[id_load,:] spzeros(nl,ng)] ./ dl
        J[id_gen, :] = [spzeros(ng,Nbus) sparse(1:ng, 1:ng, ones(ng))]
        J[Nbus+1:end, :] = [J0[id_gen,:] -sparse(1:ng, 1:ng, dg)] ./ mg
    end

    # save the frequencies at the predefined time steps
    saved_values = SavedValues(Float64, Vector{Float64})
    cb = SavingCallback((u, t, integrator) -> integrator(t, Val{1},
        idxs = 1:Nbus), saved_values, saveat = tt)
        
    # solve the swing equations
    func = ODEFunction(swing!, jac=jacobian!)
    prob = ODEProblem(func, u0, tspan, tstops = tt, callback=cb)
    solve(prob, TRBDF2())
    omega = [saved_values.saveval[i][j] for i=1:length(saved_values.saveval), j=1:Nbus]

    return saved_values.t, omega
end


function compute_phases(
    b::Vector{Float64},
    line_ids::Matrix{Int64},
    pref::Vector{Float64};
    theta::Vector{Float64} = Float64[],
    id_slack::Int64 = 1,
    tol::Float64 = 1E-6,
    maxiter::Int64 = 14,
)
    nb = maximum(line_ids)
    theta = isempty(theta) ? zeros(nb) : theta
    
    id = setdiff(1:nb, id_slack)
    nl = size(line_ids, 1)
    
    inc = sparse([line_ids[:,1]; line_ids[:,2]], [1:nl; 1:nl],
        [-ones(nl); ones(nl)])

    error = 2 * tol
    iter = 0
    while (error > tol && iter < maxiter)
        pe = -inc * (b .* sin.(inc' * theta))
        dp = pe - pref 
        v = b .* cos.(inc' * theta)
        J = inc * (v .* inc')
        x = J[id,id] \ dp[id]
        theta[id] += x
        error = maximum(abs.(dp))
        iter += 1
    end
    iter == maxiter ? println("Max iteration reached, error: ", error) : nothing
    return theta, iter
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
    idprod = findall((dm.gen .> 0.0))
    idavail = findall(dm.max_gen[idprod] .> dP) # find id large gens in the prod list
    # println(size(idprod))
    # println(idavail)
    #idprod = findall((dm.gen .> 0.0))
    id = Int64.(zeros(size(coord, 1)))  # index in the full gen list
    id2 = Int64.(zeros(size(coord, 1))) # index in the producing gen list
    for i in eachindex(coord[:, 1])
        temp = dm.idgen[idprod[idavail]]
        id[i] = idprod[idavail[argmin((dm.coord[temp, 1] .- coord[i, 1]) .^ 2 +
                                      (dm.coord[temp, 2] .- coord[i, 2]) .^ 2)]]
        id2[i] = idavail[argmin((dm.coord[temp, 1] .- coord[i, 1]) .^ 2 +
                                (dm.coord[temp, 2] .- coord[i, 2]) .^ 2)]
    end
    return id, id2
end


function load_discrete_model(
    dataname::String,
    scaling_factor::Float64
)
    data = h5read(dataname, "/")
    coord = albers_projection( data["bus_coord"] ./ (180 / pi) )
    coord ./= scaling_factor
    dm = DiscModel(
        vec(data["gen_inertia"]),
        vec(data["gen_prim_ctrl"]),
        Int64.(vec(data["gen"][:, 1])),
        findall(vec(data["bus"][:, 2]) .== 3)[1],
        coord,
        vec(data["load_freq_coef"]),
        Int64.(data["branch"][:, 1:2]),
        -1.0 ./ data["branch"][:, 4],
        vec(data["bus"][:, 3]) / 100.0,
        vec(data["bus"][:, 9]) / 180.0 * pi,
        vec(data["gen"][:, 2]) / 100.0,
        vec(data["gen"][:, 9]) / 100.0,
        size(data["bus"], 1),
        size(data["gen"], 1),
        size(data["branch"], 1),
        )
    return dm
end

