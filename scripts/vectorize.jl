function vectorize(
    isinside::BitMatrix,
    n::Array{Float64, 2},
    bx::Array{Float64, 2},
    by::Array{Float64, 2},
    p::Array{Float64, 2},
    m::Array{Float64, 2},
    d::Array{Float64, 2},
)

    # flatten the m, d, p, bx and by matrices
    # is equivalent to x = vec(reshape(x, Nx*Ny, 1))
    isinsideflat = vec(isinside)
    bxflat = vec(bx)
    byflat = vec(by)
    pflat = vec(p)
    mflat = vec(m)
    dflat = vec(d)


    # ?? neighbours along x ??
    xneigh = sparse([], [], Float64.([]), Ny * Nx, Ny * Nx)
    # ?? neighbours along x ??
    yneigh = sparse([], [], Float64.([]), Ny * Nx, Ny * Nx)
    # ?? b stuff ??
    bflat = sparse([], [], Float64.([]), Ny * Nx, Ny * Nx)
    
    #Threads.@threads for i in 1:Nx*Ny
    for i in 1:Nx*Ny
        if(isinsideflat[i])
            xneigh[i, i] = 1.0
            xneigh[i, i - Ny] = 1.0
            yneigh[i, i] = 1.0
            yneigh[i, i - 1] = 1.0
            bflat[i, i-Ny] = bxflat[i-Ny]
            bflat[i, i+Ny] = bxflat[i]
            bflat[i, i-1] = byflat[i-1]
            bflat[i, i+1] = byflat[i]
        end
    end
    
    #Threads.@threads for k in 1:size(n, 1)
    for k in 1:size(n, 1)
        # with Nx goes out of bound
        i = (Int64(n[k, 2]) - 1) * Ny + Int64(n[k, 1])
        nx = n[k, 4] 
        ny = n[k, 3]
        xneigh[i, i - Ny] = (1 + nx) # this won't store anything if fed with 0
        xneigh[i, i] = (1 - nx)
        yneigh[i, i - 1] = (1 + ny)
        yneigh[i, i] = (1 - ny)
        bflat[i, i-Ny] = (1 + nx) * bxflat[i-Ny]
        bflat[i, i+Ny] = (1 - nx) * bxflat[i]
        bflat[i, i-1] = (1 + ny) * byflat[i-1]
        bflat[i, i+1] = (1 - ny) * byflat[i]
    end  
    
    minvflat = mflat.^(-1)
    gammaflat = dflat .* minvflat
    #return isinsideflat, bxflat, byflat, pflat, mflat, dflat, bflat, xneigh, yneigh
    return isinsideflat, bxflat, byflat, pflat, minvflat, gammaflat,
        SparseMatrixCSC{Float64, Int64}(bflat),
        SparseMatrixCSC{Float64, Int64}(xneigh),
        SparseMatrixCSC{Float64, Int64}(yneigh)
end
