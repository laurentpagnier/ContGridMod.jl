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
    xi = sparse([], [], Float64.([]), Ny * Nx, Ny * Nx)

    #Threads.@threads for i in 1:Nx*Ny
    for k in 1:Nx*Ny
        if(isinsideflat[k])
            xi[k, k] = - bxflat[k] - bxflat[k-Ny] -
                byflat[k] - byflat[k-1]
            xi[k, k-Ny] = bxflat[k-Ny]
            xi[k, k+Ny] = bxflat[k]
            xi[k, k-1] = byflat[k-1]
            xi[k, k+1] = byflat[k]
        end
    end
    
    #Threads.@threads for k in 1:size(n, 1)
    for id in 1:size(n, 1)
        k = (Int64(n[id, 2]) - 1) * Ny + Int64(n[id, 1])
        ny = n[id, 3]
        nx = n[id, 4] 
        
        xi[k, k] = - (1.0 - nx) * bxflat[k] - (1.0 + nx) * bxflat[k-Ny] -
            (1.0 - ny) * byflat[k] - (1.0 + ny) * byflat[k-1]
        xi[k, k-Ny] = (1.0 + nx) * bxflat[k-Ny]
        xi[k, k+Ny] = (1.0 - nx) * bxflat[k]
        xi[k, k-1] = (1.0 + ny) * byflat[k-1]
        xi[k, k+1] = (1.0 - ny) * byflat[k]
    end  
    
    minvflat = mflat.^(-1)
    minvflat = minvflat[isgridflat]
    gammaflat = dflat[isgridflat] .* minvflat
    xi = SparseMatrixCSC{Float64, Int64}(xi[isgridflat,isgridflat])
    return isinsideflat, pflat, minvflat, gammaflat, xi
end

function vectorize_ref(
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
    xi = sparse([], [], Float64.([]), Ny * Nx, Ny * Nx)

    #Threads.@threads for i in 1:Nx*Ny
    for k in 1:Nx*Ny
        if(isinsideflat[k])
            xi[k, k] = - bxflat[k] - bxflat[k-Ny] -
                byflat[k] - byflat[k-1]
            xi[k, k-Ny] = bxflat[k-Ny]
            xi[k, k+Ny] = bxflat[k]
            xi[k, k-1] = byflat[k-1]
            xi[k, k+1] = byflat[k]
        end
    end
    
    #Threads.@threads for k in 1:size(n, 1)
    for id in 1:size(n, 1)
        k = (Int64(n[id, 2]) - 1) * Ny + Int64(n[id, 1])
        ny = n[id, 3]
        nx = n[id, 4] 
        
        xi[k, k] = - (1.0 - nx) * bxflat[k] - (1.0 + nx) * bxflat[k-Ny] -
            (1.0 - ny) * byflat[k] - (1.0 + ny) * byflat[k-1]
        xi[k, k-Ny] = (1.0 + nx) * bxflat[k-Ny]
        xi[k, k+Ny] = (1.0 - nx) * bxflat[k]
        xi[k, k-1] = (1.0 + ny) * byflat[k-1]
        xi[k, k+1] = (1.0 - ny) * byflat[k]
    end  
    
    minvflat = mflat.^(-1)
    minvflat = minvflat[isgridflat]
    gammaflat = dflat[isgridflat] .* minvflat
    xi = SparseMatrixCSC{Float64, Int64}(xi[isgridflat,isgridflat])
    return isinsideflat, pflat, minvflat, gammaflat, xi
end


