export solve_dynamics

using SparseArrays
using LinearAlgebra
using DifferentialEquations

function solve_dynamics(
    contmod::ContGridMod.ContModel,
    tspan::Tuple{Number, Number}; 
    dt::Float64 = 1e-2, # Distance of time steps at which the solution is returned
    alg = TRBDF2(),  # The solver that is passed to the solve function
    solve_kwargs::Dict = Dict(),  # Keyword arguments to the solve function
)

    # get the stable solution
    compute_stable_sol!(contmod)

    # preparations for dynamical simulation
    Nnode = contmod.mesh.Nnode
    nline = contmod.mesh.Nedge
    temp = to_mat(contmod.mesh.edge_list)
    inc = sparse([temp[:, 1]; temp[:, 2]], [1:nline; 1:nline], [-ones(nline); ones(nline)])
    u0 = [contmod.th; zeros(Nnode)]
    
    function lin_dyn!(du, u, param, t)
        dpe = contmod.p + contmod.dp - inc * (contmod.b .* (inc' * u[1:Nnode]))
        du[1:Nnode] = u[Nnode+1:end]
        du[Nnode+1:end] = (dpe - contmod.d .* u[Nnode+1:end]) ./ contmod.m
    end

    function jacobian!(J, u, param, t)
        J0 = inc * (-contmod.b .* inc')
        J[:, :] =  [spzeros(Nnode,Nnode) sparse(1:Nnode, 1:Nnode, ones(Nnode));
            J0 ./ contmod.m  sparse(1:Nnode, 1:Nnode, -contmod.d./contmod.m)]
    end

    jac_proto = spzeros(2*Nnode, 2*Nnode)
    jacobian!(jac_proto, ones(2*Nnode), [], 0.0)
    for i=1:2*Nnode
        jac_proto[i, i] += 1
    end

    # save the frequencies at the predefined time steps
    saved_values = SavedValues(Float64, Vector{Float64})
    tt = tspan[1]:dt:tspan[2]
    cb = SavingCallback((u, t, integrator) -> u[Nnode+1:end], saved_values, saveat = tt)
        
    # integrate dynamics
    func = ODEFunction(lin_dyn!, jac = jacobian!, jac_prototype = jac_proto)
    prob = ODEProblem(func, u0, tspan, tstops = tt, callback=cb)
    solve(prob, alg)
    omega = [saved_values.saveval[i][j] for i=1:length(saved_values.saveval), j=1:Nnode]

    return saved_values.t, omega
end


function solve_dynamics(
    dm::ContGridMod.DiscModel,
    tspan::Tuple{Real, Real},
    delta_p::Union{Float64, Vector{Float64}};
    faultid::Int64 = 0,
    coords::Matrix{Float64} = [NaN NaN],
    dt::Float64 = 1e-2,  # Distance of time steps at which the solution is returned
    scale_factor::Float64 = 0.0,
    tol::Float64 = 1e-10,  # Target tolerance for the Newtoon-Raphson solver
    maxiter::Int64 = 30,  # Maximum iteration for the Newton-Raphson solver
    alg=TRBDF2(),  # The solver that is passed to the solve function
    solve_kwargs::Dict = Dict(),  # Keyword arguments to the solve function
    dmin = 0.0001,
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
    temp = to_mat(dm.line_list)
    inc = sparse([temp[:, 1]; temp[:, 2]], [1:nline; 1:nline], [-ones(nline); ones(nline)])

    # create arrays of dynamical parameters
    mg = dm.m_gen[is_producing]
    dg = dm.d_gen[is_producing] + dm.d_load[id_gen]
    dl = max.(dm.d_load[id_load], dmin)

    # get the stable solution
    theta, _ = ContGridMod.compute_phases(dm.b, inc, p, theta = dm.th,
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

    jac_proto = spzeros(Nbus + ng, Nbus + ng)
    jacobian!(jac_proto, ones(Nbus + ng), 0, 0)
    for i=1:Nbus + ng
        jac_proto[i, i] += 1
    end
    
    # save the frequencies at the predefined time steps
    tt = tspan[1]:dt:tspan[2]
    saved_values = SavedValues(Float64, Vector{Float64})
    cb = SavingCallback((u, t, integrator) -> integrator(t, Val{1},
        idxs = 1:Nbus), saved_values, saveat = tt)
        
    # solve the swing equations
    func = ODEFunction(swing!, jac=jacobian!, jac_prototype = jac_proto)
    prob = ODEProblem(func, u0, tspan, tstops = tt, callback=cb)
    solve(prob, TRBDF2())
    omega = [saved_values.saveval[i][j] for i=1:length(saved_values.saveval), j=1:Nbus]

    return saved_values.t, omega
end


function compute_phases(
    b::Vector{Float64},
    inc::SparseMatrixCSC{Float64, Int64},
    pref::Vector{Float64};
    theta::Vector{Float64} = Float64[],
    id_slack::Int64 = 1,
    tol::Float64 = 1E-6,
    maxiter::Int64 = 14,
)
    theta = isempty(theta) ? zeros(nb) : theta
    nl, nb = size(inc, 2), size(inc, 1)
    id = setdiff(1:nb, id_slack)
    
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
