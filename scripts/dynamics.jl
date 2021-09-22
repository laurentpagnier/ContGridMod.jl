using SparseArrays
using LinearAlgebra


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
    isgridflat::BitArray,
    xi::SparseMatrixCSC{Float64, Int64},
    pflat::Array{Float64, 1},
    minvflat::Array{Float64, 1},
    gammaflat::Array{Float64, 1},
    th0::Array{Float64, 1};
    interval::Int64 = 10,
    Ndt::Int64 = 1000,
    dt::Float64 = 0.05
)
    println("Total time: ", dt * Ndt)
    N = sum(isgridflat)
    M = 1 + Int64(ceil(Ndt/interval))
    omegas = zeros(N, M)
    thetas = zeros(N, M)
    ts = zeros(M)
    
    x = [zeros(N); th0[isgridflat]]
    omegas[:, 1] = zeros(N)
    thetas[:, 1] = copy(th0[isgridflat])  

    ts = zeros(1 + Int64(ceil(Ndt/interval)))

    I = sparse(1:N, 1:N, ones(N))
    A = [I -dt / 2 * I;
        - dt / 2 / dx^2 * sparse(1:N, 1:N, minvflat) * xi (I + dt/2 * sparse(1:N, 1:N, gammaflat))]
    B = [I dt / 2 * I;
         dt / 2 / dx^2 * sparse(1:N, 1:N, minvflat) * xi (I - dt/2 * sparse(1:N, 1:N, gammaflat))]
    C = [zeros(N); dt * sparse(1:N, 1:N, minvflat) * pflat]

    @time begin
        for t in 1:Ndt
            x = A \ (B * x + C)
            if(mod(t,interval) == 0)
                thetas[:,Int64(t/interval) + 1] = x[1:N]
                omegas[:,Int64(t/interval) + 1] = x[N+1:end]
                ts[Int64(t/interval) + 1] = t*dt
                println("NIter: ", t, " Avg. Omega: ", sum(omegas[:, Int64(t/interval) + 1])/sum(isgridflat))
            end
        end
    end
    return ts, thetas, omegas
end


function perform_dyn_sim_backward_euler(
    isgridflat::BitArray,
    xi::SparseMatrixCSC{Float64, Int64},
    pflat::Array{Float64, 1},
    minvflat::Array{Float64, 1},
    gammaflat::Array{Float64, 1},
    th0::Array{Float64, 1};
    interval::Int64 = 10,
    Ndt::Int64 = 1000,
    dt::Float64 = 0.05
)
    println("Total time: ", dt * Ndt)
    N = sum(isgridflat)
    M = 1 + Int64(ceil(Ndt/interval))
    omegas = zeros(N, M)
    thetas = zeros(N, M)
    ts = zeros(M)
    
    x = [zeros(N); th0[isgridflat]]
    omegas[:,1] = zeros(N)
    thetas[:,1] = copy(th0[isgridflat])  

    ts = zeros(1 + Int64(ceil(Ndt/interval)))

    I = sparse(1:N, 1:N, ones(N))
    A = [I -dt * I;
        - dt / dx^2 * sparse(1:N, 1:N, minvflat) * xi (I + dt * sparse(1:N, 1:N, gammaflat))]
    B = [zeros(N); dt * sparse(1:N, 1:N, minvflat) * pflat]

    @time begin
        for t in 1:Ndt
            x = A \ (x + B)
            if(mod(t,interval) == 0)
                thetas[:,Int64(t/interval) + 1] = x[1:N]
                omegas[:,Int64(t/interval) + 1] = x[N+1:end]
                ts[Int64(t/interval) + 1] = t*dt
                println("NIter: ", t, " Avg. Omega: ", sum(omegas[:, Int64(t/interval) + 1])/sum(isgridflat))
            end
        end
    end
    return ts, thetas, omegas
end


function vectorize(
    isinside::BitMatrix,
    isborder::BitMatrix,
    n::Array{Float64, 2},
    bx::Array{Float64, 2},
    by::Array{Float64, 2},
    p::Array{Float64, 2},
    m::Array{Float64, 2},
    d::Array{Float64, 2},
)
    # flatten the m, d, p, bx and by matrices
    isinsideflat = vec(isinside)
    isgridflat = vec(isinside .| isborder)
    bxflat = vec(bx)
    byflat = vec(by)
    pflat = vec(p)
    mflat = vec(m)
    dflat = vec(d)
    Ny = size(bx,1)
    Nx = size(bx,2)

    id1 = Array{Int64,1}()
    id2 = Array{Int64,1}()
    v = Array{Float64,1}()
    for k in 1:Nx*Ny
        if(isinsideflat[k])
            append!(id1, k)
            append!(id2, k)
            append!(v, - (bxflat[k] +
                bxflat[k-Ny] +
                byflat[k] +
                byflat[k-1])
                )
            append!(id1, k)
            append!(id2, k-Ny)
            append!(v, bxflat[k-Ny])
            append!(id1, k)
            append!(id2, k+Ny)
            append!(v, bxflat[k])
            append!(id1, k)
            append!(id2, k-1)
            append!(v, byflat[k-1])
            append!(id1, k)
            append!(id2, k+1)
            append!(v, byflat[k])            
        end
    end
    
    for id in 1:size(n, 1)
        k = (Int64(n[id, 2]) - 1) * Ny + Int64(n[id, 1])
        ny = n[id, 3]
        nx = n[id, 4] 
        append!(id1, k)
        append!(id2, k)
        append!(v, - ((1.0 - nx) * bxflat[k] +
            (1.0 + nx) * bxflat[k-Ny] +
            (1.0 - ny) * byflat[k] +
            (1.0 + ny) * byflat[k-1]))
        append!(id1, k)
        append!(id2, k-Ny)
        append!(v, (1.0 + nx) * bxflat[k-Ny])
        append!(id1, k)
        append!(id2, k+Ny)
        append!(v, (1.0 - nx) * bxflat[k])
        append!(id1, k)
        append!(id2, k-1)
        append!(v, (1.0 + ny) * byflat[k-1])   
        append!(id1, k)
        append!(id2, k+1)
        append!(v, (1.0 - ny) * byflat[k])
    end
    
    xi = sparse(id1, id2, v, Ny * Nx, Ny * Nx)
    minvflat = mflat.^(-1)
    minvflat = minvflat[isgridflat]
    gammaflat = dflat[isgridflat] .* minvflat
    xi = SparseMatrixCSC{Float64, Int64}(xi[isgridflat,isgridflat])
    return isinsideflat, pflat, minvflat, gammaflat, xi
end
