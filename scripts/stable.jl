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
                #ths[:, :, Int64(t/interval)] = th
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




#=
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
        etamx = 1 - nx/2-nx^2/2
        etapx = 1 + nx/2-nx^2/2
        etamy = 1 - ny/2-ny^2/2
        etapy = 1 + ny/2-ny^2/2
        append!(id1, k)
        append!(id2, k)
        append!(v, - (etamx * bxflat[k] +
            etapx * bxflat[k-Ny] +
            etamy * byflat[k] +
            etapy * byflat[k-1]))
        append!(id1, k)
        append!(id2, k-Ny)
        append!(v, etapx * bxflat[k-Ny])
        append!(id1, k)
        append!(id2, k+Ny)
        append!(v, etamx * bxflat[k])
        append!(id1, k)
        append!(id2, k-1)
        append!(v, etapy * byflat[k-1])   
        append!(id1, k)
        append!(id2, k+1)
        append!(v, etamy * byflat[k])
    end
    
    
    xi = sparse(id1, id2, v, Ny * Nx, Ny * Nx)
    minvflat = mflat.^(-1)
    minvflat = minvflat[isgridflat]
    gammaflat = dflat[isgridflat] .* minvflat
    xi = SparseMatrixCSC{Float64, Int64}(xi[isgridflat,isgridflat])
    return isinsideflat, pflat, minvflat, gammaflat, xi
end
=#
