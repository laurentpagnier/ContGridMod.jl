using LinearAlgebra
using SparseArrays
using IterativeSolvers
include("params.jl")
include("mesh.jl")

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


function find_gen(
    dm::DiscModel,
    gps_coord::Array{Float64, 2},
    dP::Float64;
    scale_factor::Float64 = 1.0
)
    #find the the nearest generator that can "withstand" a dP fault
    coord = alberts_projection(gps_coord[:,[2;1]] ./ (180 / pi) )
    coord = coord[:,[2,1]] / scale_factor
    idprod = findall((dm.gen .> 0.0))
    idavail = findall(dm.max_gen[idprod] .> dP) # find id large gens in the prod list
    println(size(idprod))
    println(idavail)
    #idprod = findall((dm.gen .> 0.0))
    id = Int64.(zeros(size(coord,1)))  # index in the full gen list
    id2 = Int64.(zeros(size(coord,1))) # index in the producing gen list
    for i in 1:size(coord,1)
        temp = dm.idgen[idprod[idavail]]
        id[i] = idprod[idavail[argmin((dm.coord[temp,1] .- coord[i,1]).^2 +
            (dm.coord[temp,2] .- coord[i,2]).^2)]]
        id2[i] = idavail[argmin((dm.coord[temp,1] .- coord[i,1]).^2 +
            (dm.coord[temp,2] .- coord[i,2]).^2)]
    end
    return id, id2
end


function find_node(
    cm::ContModel,
    gps_coord::Array{Float64, 2};
    scale_factor::Float64 = 1.0
)
    coord = alberts_projection(gps_coord[:,[2;1]] ./ (180 / pi) )
    coord = coord[:,[2,1]] / scale_factor

    id = Int64.(zeros(size(coord,1))) # index in the full gen list
    for i in 1:size(coord,1)
        id[i] = argmin((cm.coord[cm.isgrid,1] .- coord[i,1]).^2 +
            (cm.coord[cm.isgrid,2] .- coord[i,2]).^2)
    end
    return id
end

function find_node2(
    cm::ContModel,
    coord::Array{Float64, 2}
)
    id = Int64.(zeros(size(coord,1))) # index in the full gen list
    for i in 1:size(coord,1)
        id[i] = argmin((cm.coord[cm.isgrid,1] .- coord[i,1]).^2 +
            (cm.coord[cm.isgrid,2] .- coord[i,2]).^2)
    end
    return id
end
