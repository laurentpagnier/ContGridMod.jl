function perform_dyn_sim(
    isinside::BitMatrix,
    n::Array{Float64, 2},
    bx::Array{Float64, 2},
    by::Array{Float64, 2},
    p::Array{Float64, 2},
    m::Array{Float64, 2},
    d::Array{Float64, 2},
    th0::Array{Float64, 2};
    interval::Int64 = 100,
    Ndt::Int64 = 15000,
    dt = 0.0001
)

    temp = findall(isinside)
    idin = Int64.(zeros(size(temp, 1), 2))
    for i in 1:size(temp,1)
        idin[i,:] = [temp[i][1], temp[i][2]]
    end
    
    gamma = d ./ m

    omegas = zeros(Ny, Nx, 1 + Int64(ceil(Ndt/interval)))
    thetas = zeros(Ny, Nx, 1 + Int64(ceil(Ndt/interval)))
    th = copy(th0)
    th_new = zeros(Ny, Nx)
    th_old = copy(th0)
    omegas[:,:,1] = (th - th_old) / dt
    thetas[:,:,1] = th

    ts = zeros(1 + Int64(ceil(Ndt/interval)))
    chi = 1.0 .+ gamma*dt ./ 2.0

    b = zeros(Ny, Nx)
    Threads.@threads for k in 1:size(idin, 1)
        i = idin[k, 1]
        j = idin[k, 2]
        b[i, j] = (by[i-1, j] + by[i, j] + bx[i, j-1] + bx[i, j])
    end
            
    Threads.@threads for k in 1:size(n, 1)
        i = Int64(n[k,1])
        j = Int64(n[k,2])
        nx = n[k,4]
        ny = n[k,3]
        b[i, j] = (1.0 + ny) * by[i-1, j] +
            (1.0 - ny) * by[i, j] +
            (1.0 + nx) * bx[i, j-1] +
            (1.0 - nx) * bx[i, j]
    end 

    @time begin
        for t in 1:Ndt
            
            Threads.@threads for k in 1:size(idin, 1)
                i = idin[k, 1]
                j = idin[k, 2]
                th_new[i, j] = (2.0 - b[i, j] / m[i, j] * dt^2 / dx^2) / chi[i, j] * th[i, j] + 
                    ( gamma[i, j] * dt / 2.0 - 1.0 ) / chi[i, j] * th_old[i, j] + 
                    dt^2 / dx^2 / chi[i, j] / m[i, j] * (
                    by[i, j] * th[i+1, j] +
                    by[i-1, j] * th[i-1, j] +
                    bx[i, j] * th[i, j+1] +
                    bx[i, j-1] * th[i, j-1]
                    ) + dt^2 / chi[i, j] / m[i, j] * p[i, j]
            end
            
            Threads.@threads for k in 1:size(n, 1)
                i = Int64(n[k, 1])
                j = Int64(n[k, 2])
                nx = n[k, 4]
                ny = n[k, 3]
                th_new[i, j] = (2.0 - b[i, j] / m[i, j] * dt^2 / dx^2) / chi[i, j] * th[i, j] + 
                    ( gamma[i, j] * dt / 2.0 - 1.0 ) / chi[i, j] * th_old[i, j] + 
                    dt^2 / dx^2 / chi[i, j] / m[i, j] * (
                    (1 - ny) * by[i, j] * th[i+1, j] +
                    (1 + ny) * by[i-1, j] * th[i-1, j] +
                    (1 - nx) * bx[i, j] * th[i, j+1] +
                    (1 + nx) * bx[i, j-1] * th[i, j-1]
                    ) + dt^2 / chi[i, j] / m[i, j] * p[i, j]
            end  

            if(mod(t, interval) == 0)
                omegas[:,:,Int64(t/interval) + 1] = (th_new - th) / dt
                thetas[:,:,Int64(t/interval) + 1] = th_new
                ts[Int64(t/interval) + 1] = t * dt
                println("NIter: ", t, " Omega: ", sum(omegas[isinside,Int64(t/interval) + 1]) / sum(isinside))
            end
            
            th_old = copy(th)
            th = copy(th_new)
        end
    end
    return ts, thetas, omegas
end



function perform_dyn_sim_vec(
    isgridflat::BitArray,
    xi::SparseMatrixCSC{Float64, Int64},
    pflat::Array{Float64, 1},
    minvflat::Array{Float64, 1},
    gammaflat::Array{Float64, 1},
    th0::Array{Float64, 1};
    interval::Int64 = 100,
    Ndt::Int64 = 15000,
    dt = 0.0001
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



function perform_dyn_sim_vec_crank_nicolson(
    isgridflat::BitArray,
    xi::SparseMatrixCSC{Float64, Int64},
    pflat::Array{Float64, 1},
    minvflat::Array{Float64, 1},
    gammaflat::Array{Float64, 1},
    th0::Array{Float64, 1};
    interval::Int64 = 10,
    Ndt::Int64 = 1000,
    dt = 0.05
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



function perform_dyn_sim_vec_backward_euler(
    isgridflat::BitArray,
    xi::SparseMatrixCSC{Float64, Int64},
    pflat::Array{Float64, 1},
    minvflat::Array{Float64, 1},
    gammaflat::Array{Float64, 1},
    th0::Array{Float64, 1};
    interval::Int64 = 10,
    Ndt::Int64 = 1000,
    dt = 0.05
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
