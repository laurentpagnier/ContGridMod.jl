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
    omegas[:,:,1] = (th - th_old) / dt;
    thetas[:,:,1] = th;
    #println(sum(p))
    ts = zeros(1 + Int64(ceil(Ndt/interval)))
    chi = 1.0 .+ gamma*dt ./ 2.0
    #println(chi)
    @time begin
        for t in 1:Ndt
            
            Threads.@threads for k in 1:size(idin, 1)
                i = idin[k, 1]
                j = idin[k, 2]
                bij = ( by[i-1, j] + by[i, j] + bx[i, j-1] + bx[i, j] )
                th_new[i, j] = (2.0 - bij / m[i, j] * dt^2 / dx^2) / chi[i, j] * th[i, j] + 
                    ( gamma[i, j] * dt / 2.0 - 1.0 ) / chi[i, j] * th_old[i, j] + 
                    dt^2 / dx^2 / chi[i, j] / m[i, j] * (
                    by[i, j] * th[i+1, j] +
                    by[i-1, j] * th[i-1, j] +
                    bx[i, j] * th[i, j+1] +
                    bx[i, j-1] * th[i, j-1]
                    ) + dt^2 / chi[i, j] / m[i, j] * p[i, j]
            end
            
            Threads.@threads for k in 1:size(n,1)
                i = Int64(n[k,1])
                j = Int64(n[k,2])
                nx = n[k,4]
                ny = n[k,3]
                bij = (1.0 + ny) * by[i-1, j] + (1.0 - ny) * by[i, j] +
                    (1.0 + nx) * bx[i, j-1] + (1.0 - nx) * bx[i, j]
                th_new[i, j] = (2.0 - bij / m[i, j] * dt^2 / dx^2) / chi[i, j] * th[i, j] + 
                    ( gamma[i, j] * dt / 2.0 - 1.0 ) / chi[i, j] * th_old[i, j] + 
                    dt^2 / dx^2 / chi[i, j] / m[i, j] * (
                    (1 - ny) * by[i, j] * th[i+1, j] +
                    (1 + ny) * by[i-1, j] * th[i-1, j] +
                    (1 - nx) * bx[i, j] * th[i, j+1] +
                    (1 + nx) * bx[i, j-1] * th[i, j-1]
                    ) + dt^2 / chi[i, j] / m[i, j] * p[i, j]
            end  

            if(mod(t, interval) == 0)
                println("NIter: ", t)
                omegas[:,:,Int64(t/interval) + 1] = (th_new - th) / dt
                thetas[:,:,Int64(t/interval) + 1] = th_new
                ts[Int64(t/interval) + 1] = t * dt
                
                Threads.@threads for k in 1:size(n,1)
                    i = Int64(n[k,1])
                    j = Int64(n[k,2])
                    nx = n[k,4]
                    ny = n[k,3]
                    if(ny == -1)
                        #println([th[i-1, j] th[i, j] th[i+1, j] th_new[i-1, j] th_new[i, j] th_new[i+1, j] ])
                    end
                end  
                
            end
            
            th_old = copy(th)
            th = copy(th_new)
        end
    end
    return ts, thetas, omegas
end
