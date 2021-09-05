function compute_stable_sol(
    isinside::BitMatrix,
    n::Array{Float64, 2},
    bx::Array{Float64, 2},
    by::Array{Float64, 2},
    p::Array{Float64, 2};
    interval::Int64 = 1000,
    Niter::Int64 = 15000,
    tol::Float64 = 1e-7
)

    temp = findall(isinside)
    idin = Int64.(zeros(size(temp, 1), 2))
    for i in 1:size(temp,1)
        idin[i,:] = [temp[i][1], temp[i][2]]
    end

    th_old = zeros(Ny, Nx)
    th = zeros(Ny, Nx)
    
    @time begin
        for t in 1:Niter
            if(mod(t, interval) == 0)
                temp = copy(th)
            end

            Threads.@threads for k in 1:size(idin, 1)
                i = idin[k, 1]
                j = idin[k, 2]
                bij = (by[i-1,j] + by[i,j] + bx[i-1, j] + bx[i,j])
                th[i,j] = (
                    by[i,j] * th_old[i+1,j] +
                    by[i-1, j] * th_old[i-1, j] + 
                    bx[i,j] * th_old[i,j+1] +
                    bx[i,j-1] * th_old[i,j-1] +
                    dx^2 * p[i,j]
                    ) / bij
            end
            
            Threads.@threads for k in 1:size(n,1)
                i = Int64(n[k,1])
                j = Int64(n[k,2])
                
                nx = n[k,4]
                ny = n[k,3]
                bij = (1 + ny) * by[i-1, j] + (1 - ny) * by[i, j] +
                    (1 + nx) * bx[i, j-1] + (1 - nx) * bx[i, j]

                th[i, j] = (
                    (1 + ny) * by[i-1, j] * th_old[i-1, j] +
                    (1 - ny) * by[i, j] * th_old[i+1, j] + 
                    (1 + nx) * bx[i, j-1] * th_old[i, j-1] +
                    (1 - nx) * bx[i, j] * th_old[i,j+1] +
                    dx^2 * p[i,j]
                    ) / bij
            end
            if(mod(t, interval) == 0)
                println( [t maximum(abs.(th - temp))] )
                if( maximum(abs.(th - temp)) < tol )
                    break
                end
            end
            th_old = copy(th)
        end
    end
    return th
end
