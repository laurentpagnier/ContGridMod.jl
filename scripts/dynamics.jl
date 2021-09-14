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
    isinsideflat::BitArray,
    bxflat::Array{Float64, 1},
    byflat::Array{Float64, 1},
    xneigh::SparseMatrixCSC{Float64, Int64},
    yneigh::SparseMatrixCSC{Float64, Int64},
    bflat::SparseMatrixCSC{Float64, Int64},
    pflat::Array{Float64, 1},
    minvflat::Array{Float64, 1},
    gammaflat::Array{Float64, 1},
    th0::Array{Float64, 1};
    interval::Int64 = 100,
    Ndt::Int64 = 15000,
    dt = 0.0001
)
    println("Total time: ", dt * Ndt)

    omegas = zeros(Ny * Nx, 1 + Int64(ceil(Ndt/interval)))
    thetas = zeros(Ny * Nx, 1 + Int64(ceil(Ndt/interval)))
    th_new = zeros(Ny * Nx)
    th_old = copy(th0)
    th = copy(th0)
    omegas[:,1] = zeros(size(th))
    thetas[:,1] = copy(th0)  

    ts = zeros(1 + Int64(ceil(Ndt/interval)))
    chi = 1 ./ (1 .+ gammaflat*dt/2)

    @time begin
        for t in 1:Ndt
            th_new = 2 * chi .* th - dt^2/dx^2 * minvflat .* chi .* (xneigh * bxflat + yneigh * byflat) .* th -
            (1 .- gammaflat * dt / 2) .* chi .* th_old + (dt^2 / dx^2) * chi .* minvflat .* (bflat * th) +
            dt^2 * chi .* minvflat .* pflat
            if(mod(t,interval) == 0)
                omegas[:,Int64(t/interval) + 1] = (th_new-th) / dt
                thetas[:,Int64(t/interval) + 1] = th_new
                ts[Int64(t/interval) + 1] = t * dt
                print("NIter: ", t, " Avg. Omega: ", mean(omegas[isinsideflat, Int64(t/interval) + 1]), "\n")
            end
            th_old = copy(th)
            th = copy(th_new)
        end
    end
    # Rewrite omegas and thetas in 2d format
    omegasre = zeros(Ny,Nx,1 + Int64(ceil(Ndt/interval)))
    thetasre = zeros(Ny,Nx,1 + Int64(ceil(Ndt/interval)))
    for i=1:1 + Int64(ceil(Ndt/interval))
        for j=1:Nx*Ny
            omegasre[(j-1) % Ny + 1, (j-1) ÷ Ny + 1, i] = omegas[j, i]
            thetasre[(j-1) % Ny + 1, (j-1) ÷ Ny + 1, i] = omegas[j, i]
        end
    end
    return ts, thetasre, omegasre
end

