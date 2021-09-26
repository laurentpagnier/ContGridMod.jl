using SparseArrays
using LinearAlgebra
using IterativeSolvers

function perform_dyn_sim(
    isgridflat::BitArray,
    xi::SparseMatrixCSC{Float64, Int64},
    pflat::Array{Float64, 1},
    minvflat::Array{Float64, 1},
    gammaflat::Array{Float64, 1},
    th0::Array{Float64, 1};
    interval::Int64,
    Ndt::Int64,
    dt::Float64,
    method::String="crank-nicolson"
)

    if(method == "crank-nicolson")
        return ts, thetas, omegas = perform_dyn_sim_crank_nicolson(isgridflat, xi,
            pflat, minvflat, gammaflat, th0; interval, Ndt, dt)  
    elseif(method == "backward-euler")
        return ts, thetas, omegas = perform_dyn_sim_backward_euler(isgridflat, xi,
            pflat, minvflat, gammaflat, th0; interval, Ndt, dt)    
    elseif(method == "forward")
        return ts, thetas, omegas = perform_dyn_sim_forward(isgridflat, xi,
            pflat, minvflat, gammaflat, th0; interval, Ndt, dt)
    else
        println("Method not found, '", method, "' is not implemented.")
    end

end


function perform_dyn_sim_forward(
    isgridflat::BitArray,
    xi::SparseMatrixCSC{Float64, Int64},
    pflat::Array{Float64, 1},
    minvflat::Array{Float64, 1},
    gammaflat::Array{Float64, 1},
    th0::Array{Float64, 1};
    interval::Int64 = 100,
    Ndt::Int64 = 15000,
    dt::Float64 = 0.0001
)
    println("Total time: ", dt * Ndt)
    N = sum(isgridflat)
    M = 1 + Int64(ceil(Ndt/interval))
    omegas = zeros(N, M)
    thetas = zeros(N, M)
    th_new = zeros(N)
    ts = zeros(M)
    
    th_old = copy(th0[isgridflat])
    th = copy(th0[isgridflat])
    omegas[:,1] = zeros(size(th))
    thetas[:,1] = copy(th0[isgridflat])  

    chi = 1 ./ (1 .+ gammaflat*dt/2)

    A = (2 * chi - dt^2/dx^2 * minvflat)
    B = (1 .- gammaflat * dt / 2) .* chi
    C = dt^2 * chi .* minvflat .* pflat
    
    @time begin
        for t in 1:Ndt
            th_new = 2 * chi .* th +
                dt^2 / dx^2 .* minvflat .* chi .* (xi * th) -
                B .* th_old + C
            if(mod(t,interval) == 0)
                omegas[:,Int64(t/interval) + 1] = (th_new-th) / dt
                thetas[:,Int64(t/interval) + 1] = th_new
                ts[Int64(t/interval) + 1] = t * dt
                println("NIter: ", t, " Avg. Omega: ", sum(omegas[:, Int64(t/interval) + 1])/sum(isgridflat))
            end
            th_old = copy(th)
            th = copy(th_new)
        end
    end

    return ts, thetas, omegas
end


function perform_dyn_sim_crank_nicolson(
    isgrid::BitArray,
    xi::SparseMatrixCSC{Float64, Int64},
    p::Array{Float64, 1},
    minv::Array{Float64, 1},
    gamma::Array{Float64, 1},
    th0::Array{Float64, 1};
    interval::Int64 = 10,
    Ndt::Int64 = 1000,
    dt::Float64 = 0.05
)
    println("Total time: ", dt * Ndt)
    N = sum(isgrid)
    M = 1 + Int64(ceil(Ndt/interval))
    omegas = zeros(N, M)
    thetas = zeros(N, M)
    ts = zeros(M)
    
    x = [th0[isgrid]; zeros(N)]
    omegas[:, 1] = zeros(N)
    thetas[:, 1] = copy(th0[isgrid])  

    ts = zeros(1 + Int64(ceil(Ndt/interval)))

    I = sparse(1:N, 1:N, ones(N))
    A = [I -dt / 2 * I;
        - dt / 2 / dx^2 * sparse(1:N, 1:N, minv) * xi (I + dt/2 * sparse(1:N, 1:N, gamma))]
    B = [I dt / 2 * I;
         dt / 2 / dx^2 * sparse(1:N, 1:N, minv) * xi (I - dt/2 * sparse(1:N, 1:N, gamma))]
    C = [zeros(N); dt * sparse(1:N, 1:N, minv) * p]

    @time begin
        for t in 1:Ndt
            #x = A \ (B * x + C) # way slower when dx -> 0
            gmres!(x, A , B * x + C)
            if(mod(t,interval) == 0)
                thetas[:,Int64(t/interval) + 1] = x[1:N]
                omegas[:,Int64(t/interval) + 1] = x[N+1:end]
                ts[Int64(t/interval) + 1] = t*dt
                println("NIter: ", t, " Avg. Omega: ", sum(omegas[:, Int64(t/interval) + 1])/sum(isgrid))
            end
        end
    end
    return ts, thetas, omegas
end


function perform_dyn_sim_backward_euler(
    isgrid::BitArray,
    xi::SparseMatrixCSC{Float64, Int64},
    p::Array{Float64, 1},
    minv::Array{Float64, 1},
    gamma::Array{Float64, 1},
    th0::Array{Float64, 1};
    interval::Int64 = 10,
    Ndt::Int64 = 1000,
    dt::Float64 = 0.05
)
    println("Total time: ", dt * Ndt)
    N = sum(isgrid)
    M = 1 + Int64(ceil(Ndt/interval))
    omegas = zeros(N, M)
    thetas = zeros(N, M)
    ts = zeros(M)
    
    x = [th0[isgrid]; zeros(N)]
    omegas[:,1] = zeros(N)
    thetas[:,1] = copy(th0[isgrid])  

    ts = zeros(1 + Int64(ceil(Ndt/interval)))

    I = sparse(1:N, 1:N, ones(N))
    A = [I -dt * I;
        - dt / dx^2 * sparse(1:N, 1:N, minv) * xi (I + dt * sparse(1:N, 1:N, gamma))]
    B = [zeros(N); dt * sparse(1:N, 1:N, minv) * p]

    @time begin
        for t in 1:Ndt
            #x = A \ (x + B)
            gmres!(x, A , x + B)
            if(mod(t,interval) == 0)
                thetas[:,Int64(t/interval) + 1] = x[1:N]
                omegas[:,Int64(t/interval) + 1] = x[N+1:end]
                ts[Int64(t/interval) + 1] = t*dt
                println("NIter: ", t, " Avg. Omega: ", sum(omegas[:, Int64(t/interval) + 1])/sum(isgrid))
            end
        end
    end
    return ts, thetas, omegas
end
