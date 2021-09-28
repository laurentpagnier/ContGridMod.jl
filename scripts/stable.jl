function compute_stable_sol!(
    contmod = ContModel;
    interval::Int64 = 1000,
    Niter::Int64 = 15000,
    tol::Float64 = 1e-7
)
    id1 = Array{Int64,1}()
    id2 = Array{Int64,1}()
    v = Array{Float64,1}()
    b = zeros(size(contmod.isinside))
    N = contmod.Nx * contmod.Ny
    for k in 1:N
        if(contmod.isinside[k])
            append!(id1, k)
            append!(id2, k-contmod.Ny)
            append!(v, contmod.bx[k-contmod.Ny])
            append!(id1, k)
            append!(id2, k+contmod.Ny)
            append!(v, contmod.bx[k])
            append!(id1, k)
            append!(id2, k-1)
            append!(v, contmod.by[k-1])
            append!(id1, k)
            append!(id2, k+1)
            append!(v, contmod.by[k])  
            b[k] = (contmod.bx[k] +
                contmod.bx[k-contmod.Ny] +
                contmod.by[k] +
                contmod.by[k-1])
        end
    end
    
    for id in 1:size(contmod.n, 1)
        k = (Int64(contmod.n[id, 2]) - 1) * contmod.Ny + Int64(contmod.n[id, 1])
        ny = contmod.n[id, 3]
        nx = contmod.n[id, 4] 
        etamx = 1 - nx/2 - nx^2/2
        etapx = 1 + nx/2 - nx^2/2
        etamy = 1 - ny/2 - ny^2/2
        etapy = 1 + ny/2 - ny^2/2
        append!(id1, k)
        append!(id2, k-contmod.Ny)
        append!(v, etapx * contmod.bx[k-contmod.Ny])
        append!(id1, k)
        append!(id2, k+contmod.Ny)
        append!(v, etamx * contmod.bx[k])
        append!(id1, k)
        append!(id2, k-1)
        append!(v, etapy * contmod.by[k-1])   
        append!(id1, k)
        append!(id2, k+1)
        append!(v, etamy * contmod.by[k])
        b[k] = (etamx * contmod.bx[k] +
            etapx * contmod.bx[k-contmod.Ny] +
            etamy * contmod.by[k] +
            etapy * contmod.by[k-1])
    end
    
    xi = sparse(id1, id2, v, N, N)
    p2 = dx^2 * contmod.p
    th = zeros(sum(contmod.isgrid))
    xi = xi[contmod.isgrid, contmod.isgrid]
    b = b[contmod.isgrid]
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
    contmod.th[contmod.isgrid] = th
end
