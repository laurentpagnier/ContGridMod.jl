using LinearAlgebra
using SparseArrays
using IterativeSolvers

export find_gen, find_node, find_node2, radau5, NRsolver, disc_dynamics
#=
TODO: Discuss how the disc_dynamics should be treated, should we reorder it or expect the data to be in order?
=#

function disc_dynamics(
        dm::DiscModel,
        dP::Float64;
        faultid::Int64 = 0,
        coords::Array{Float64, 2} = [NaN NaN],
        interval::Int64 = 100,
        Ndt::Int64 = 25000,
        dt::Float64 = 5e-3,
        scale_factor::Float64 = 0.0
)
    #=
    This function allows to either specify the GPS coordinates of where the fault occurs or the ID of the generator.
    If both are given the ID will be used. If the GPS coordinates are given, the scale_factor must be given.
    =#
    if (coords == [NaN NaN] && faultid == 0)
        throw(ErrorException("Either the coordinates or the ID of the faulty generator need to be specified."))
    end
    if (scale_factor == 0.0 && faultid == 0)
        throw(ErrorException("If the coordinates are given, scale_factor must be specified."))
    end
    if (faultid == 0)
        ~, tmp = find_gen(dm, coords, dP, scale_factor=scale_factor)
        faultid = tmp[1]
    end
    # Data preperation
    Nbus = length(dm.load)
    is_producing = dm.gen .> 0
    id_gen = dm.idgen[is_producing]
    id_load = setdiff(1:Nbus, id_gen)
    ng = length(id_gen)
    nl = length(id_load)
    g = zeros(Nbus)
    g[dm.idgen] .= dm.gen
    p = -dm.load + g
    p .-= sum(p) / Nbus
    edges = Int64.(zeros(size(dm.idb)))
    line_start = dm.idb[:, 1]
    line_end = dm.idb[:, 2]
    # bus reordering: generator buses first 
    for i in 1:ng
        edges[line_start .== id_gen[i], 1] .= i
        edges[line_end .== id_gen[i], 2] .= i
    end
    for i in 1:nl
        edges[line_start .== id_load[i], 1] .= i + ng
        edges[line_end .== id_load[i], 2] .= i + ng
    end
    nline = size(dm.idb,1)
    id = edges
    inc = sparse([id[:,1]; id[:,2]], [1:nline; 1:nline], [-ones(nline);ones(nline)])

    # set initial conditions
    mg = dm.mg[is_producing]
    dg = dm.dg[is_producing]
    dg += dm.dl[id_gen]
    dl = dm.dl[id_load]
    pg = p[id_gen]
    pl = p[id_load]
    p = [pg; pl]
    # get the stable solution
    b2 = - im .* dm.bline
    Ybus = conj(inc * sparse(1:nline, 1:nline, b2) * inc')
    q = zeros(Nbus)
    V = ones(Nbus)
    theta = zeros(Nbus)

    V, theta, iter = NRsolver(Ybus, V, theta, -p, q, Array{Int64,1}([]), 1, tol = 1E-30, maxiter = 10)
    thg = theta[1:ng]
    thl = theta[ng+1:end]
    og = zeros(ng);
    ts = zeros(Int64.(floor(Ndt / interval))+1)
    omegas = zeros(ng, Int64.(floor(Ndt / interval))+1)
    k = 1
    dp = zeros(ng)
    dp[faultid] = -9.0
    for i in 1:Ndt
        y = radau5(og, thg, thl, mg, dg, dl, pg+dp, pl, inc, dm.bline, dt, maxiter = 14, tol = 1E-6)
        if(mod(i, interval) == 0)
            println(i, " / ", Ndt, " (", floor(100*i/Ndt), "%)")
                omegas[:, k+1] = y[1 : ng]
                ts[k+1] = i*dt
                k += 1
            end
        og = y[1:ng]
        thg = y[ng + 1 : 2 * ng]
        thl = y[2 * ng + 1 : end]
    end
    return ts, omegas
end

function radau5(
    og::Array{Float64, 1},
    thg::Array{Float64, 1},
    thl::Array{Float64, 1},
    mg::Array{Float64, 1},
    dg::Array{Float64, 1},
    dl::Array{Float64, 1},
    pg::Array{Float64, 1},
    pl::Array{Float64, 1},
    incidence_mat::SparseMatrixCSC{Float64, Int64},
    line_susceptance::Array{Float64, 1},
    dt::Float64 = 0.0001;
    maxiter::Int64 = 14,
    tol::Float64 = 1E-6
)

    #=
    This method is adapted from its version on the
    Pantagruel repository (https://doi.org/10.5281/zenodo.2642175).
    For information on Radau methods, see, for instance,
    E. Hairer and G. Wanner, Stiff differential equations solved by Radau
    methods, J. Comput. Appl. Math. 11(1-2): 93-111 (1999)
    =#

    nline = length(line_susceptance)
    Ap = abs.(incidence_mat)
    B = sparse(1:nline, 1:nline, line_susceptance)
    ng = length(og) # number of generators
    nb = ng + length(thl) # number of buses
    nvar = ng + nb # number of variables
    ns = 3
    
    # butcher tableau for radau 5 
    a = [(88.0 - 7.0 * sqrt(6.0)) / 360.0 (296.0 - 169.0 * sqrt(6.0)) / 1800.0 (-2.0 + 3.0 * sqrt(6.0)) / 225.0;
        (296.0 + 169.0 * sqrt(6.0)) / 1800.0  (88.0 + 7.0 * sqrt(6.0)) / 360.0 (-5.0 - 3.0 * sqrt(6.0)) / 225.0;
        (16.0 - sqrt(6.0)) / 36.0 (16.0 + sqrt(6.0)) / 36.0 1.0 / 9.0]
        
    # diagonalisation of A
    res = eigen(a)
    lambda = 1.0 ./ res.values
    Tr = res.vectors
    Tl = inv(Tr)
    
    y0 = [og;
        thg;
        thl]
    Id_var = sparse(1:nvar, 1:nvar, ones(nvar))
    M1 = kron(Tr, Id_var)
    M2 = kron(lambda .* Tl, Id_var)
    
    mat1 = incidence_mat[1:ng, :] * B
    mat2 = incidence_mat[ng+1:end,:] * B
    dtheta = incidence_mat' * y0[ng+1:end]
    
    # Jacobian matrix
    M = sparse(1:nvar, 1:nvar, [mg ; ones(ng); dl]) # Mass matrix
    J0 = Ap * sparse(1:nline, 1:nline, line_susceptance .* cos.(dtheta)) * Ap'
    J0 = J0 - sparse(1:nb, 1:nb, 2 * diag(J0))
    
    J = [sparse(1:ng, 1:ng, -dg) J0[1:ng,:];
        sparse(1:ng, 1:ng, ones(ng)) SparseMatrixCSC{Float64, Int64}(sparse([], [], [], ng, nb));
        SparseMatrixCSC{Float64, Int64}(sparse([], [], [], nb - ng, ng)) J0[ng+1:end,:]]

    P = [pg; zeros(ng); pl]
   
    Y = vec(repeat(y0, ns, 1))
    F = Complex.(zeros(length(y0) * ns))
    W = Complex.(zeros(length(y0) * ns))
    dW = Complex.(zeros(length(y0) * ns))
    iter = 0
    error = 2.0 * tol
    # Newton method
	while(error > tol && iter < maxiter)
        for s in 1:ns # for the different stages
            id = nvar * (s - 1) + 1 : nvar * s 
            F[id] = M * (Y[id] - y0)
            for r in 1:ns
                id2 = nvar * (r - 1) + ng + 1 : nvar * r
                id3 = nvar * (r - 1) + 1 : nvar * (r - 1) + ng
                dtheta = incidence_mat' * Y[id2]
                F[id] -= dt * a[s, r] * (P + [-dg .* Y[id3] -
                    mat1 * sin.(dtheta); 
                    Y[id3]; -mat2 * sin.(dtheta)]);
            end
        end
        F2 = - M2 * F
        for s in 1:ns
            id = nvar * (s - 1) + 1 : nvar * s
            dW[id] = (lambda[s] * M - dt * J) \ F2[id]
            #dW[id] = gmres((lambda[s] * M - dt * J), F2[id])
        end
        W += dW
        error = maximum( abs.(M1 * dW) )
        Y = M1 * W + repeat(y0, ns, 1)
        # can probably be replace by  something like Y += ...
        iter += 1
	end
    if(iter == maxiter)
        println("Max iteration reached, error: ", error)
    end
    return real.(Y[(ns - 1) * nvar + 1 : end]) 
end


function NRsolver(
    Ybus::SparseMatrixCSC{ComplexF64, Int64},
    V::Array{Float64, 1},
    theta::Array{Float64, 1},
    p::Array{Float64, 1},
    q::Array{Float64, 1},
    idpq::Array{Int64, 1},
    id_slack::Int64;
    tol::Float64 = 1E-6,
    maxiter::Int64 = 14
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
    error  = 2 * tol
    iter = 0  
    id = [1 : id_slack-1; id_slack + 1 : nb]
    while(error > tol && iter < maxiter)
        Vc = V .* exp.(im * theta)
        S = Vc .* conj(Ybus * Vc)
        dPQ = [real(S[id]) - p[id]; imag(S[idpq]) - q[idpq]];
        dsdth = -im * sparse(1:nb, 1:nb, Vc) * conj(Ybus) * sparse(1:nb, 1:nb, conj(Vc)) +
            im * sparse(1:nb, 1:nb, Vc .* conj(Ybus * Vc))
        dsdv = sparse(1:nb, 1:nb, Vc) * conj(Ybus) * sparse(1:nb, 1:nb, exp.(-im * theta)) +
            sparse(1:nb, 1:nb, exp.(im * theta) .* conj(Ybus * Vc))
        J = [real(dsdth[id, id]) real(dsdv[id, idpq]); imag(dsdth[idpq, id]) imag(dsdv[idpq, idpq])];
        x = J \ dPQ
        theta[id] = theta[id] - x[1:nb-1]
        if(!isempty(idpq))
            V[idpq] -= x[n:end]
        end
        error = maximum(abs.(dPQ))
        iter += 1
    end
	if(iter == maxiter)
        println("Max iteration reached, error: ", error)
    end
    return V, theta, iter 
end



