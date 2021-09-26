function compute_stable_sol(
    isinside::BitVector,
    n::Array{Float64, 2},
    bx::Array{Float64, 1},
    by::Array{Float64, 1},
    p::Array{Float64, 1};
    interval::Int64 = 1000,
    Niter::Int64 = 15000,
    tol::Float64 = 1e-7
)

    id1 = Array{Int64,1}()
    id2 = Array{Int64,1}()
    v = Array{Float64,1}()
    b = zeros(size(isinside))
    for k in 1:Nx*Ny
        if(isinside[k])
            append!(id1, k)
            append!(id2, k-Ny)
            append!(v, bx[k-Ny])
            append!(id1, k)
            append!(id2, k+Ny)
            append!(v, bx[k])
            append!(id1, k)
            append!(id2, k-1)
            append!(v, by[k-1])
            append!(id1, k)
            append!(id2, k+1)
            append!(v, by[k])  
            b[k] = (bx[k] +
                bx[k-Ny] +
                by[k] +
                by[k-1])
        end
    end
    
    for id in 1:size(n, 1)
        k = (Int64(n[id, 2]) - 1) * Ny + Int64(n[id, 1])
        ny = n[id, 3]
        nx = n[id, 4] 
        etamx = 1 - nx/2-nx^2/2
        etapx = 1 + nx/2-nx^2/2
        etamy = 1 - ny/2-ny^2/2
        etapy = 1 + ny/2-ny^2/2
        append!(id1, k)
        append!(id2, k-Ny)
        append!(v, etapx * bx[k-Ny])
        append!(id1, k)
        append!(id2, k+Ny)
        append!(v, etamx * bx[k])
        append!(id1, k)
        append!(id2, k-1)
        append!(v, etapy * by[k-1])   
        append!(id1, k)
        append!(id2, k+1)
        append!(v, etamy * by[k])
        b[k] = (etamx * bx[k] +
            etapx * bx[k-Ny] +
            etamy * by[k] +
            etapy * by[k-1])
    end
    
    xi = sparse(id1, id2, v, Ny * Nx, Ny * Nx)
    p2 = dx^2 * p
    th = zeros(sum(isgrid))
    xi = xi[isgrid, isgrid]
    b = b[isgrid]
    @time begin
        for t in 1:Niter
            if(mod(t, interval) == 0)
                temp = copy(th)
            end
            th = (xi * th + p2) ./ b            
            if(mod(t, interval) == 0)
                println( [t maximum(abs.(th - temp))] )
                if( maximum(abs.(th - temp)) < tol )
                    break
                end
            end
        end
    end
    th2 = zeros(size(isgrid))
    th2[isgrid] = th
    return th2
end
