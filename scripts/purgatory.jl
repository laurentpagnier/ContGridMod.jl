function get_grid(border::Array{Float64,2}, dx::Float64)

    xlim = [minimum(border[:,1])-2*dx maximum(border[:,1])+2*dx]
    ylim = [minimum(border[:,2])-2*dx maximum(border[:,2])+2*dx]
    xrange = collect(xlim[1]:dx:xlim[2])
    yrange = collect(ylim[1]:dx:ylim[2])

    Nx = length(xrange)
    Ny = length(yrange)

    isborder = Bool.(zeros(Ny,Nx))
    isinside = Bool.(zeros(Ny,Nx))
    isgrid = Bool.(zeros(Ny,Nx))

    for j=1:Nx
        for i=1:Ny
            p = [xrange[j] yrange[i]]
            tmp = inPolygon(p,border)[1]
            isinside[i,j] = tmp
            isgrid[i, j] = tmp
        end
    end
    # check if alone 
    for j=2:Nx-1
        for i=2:Ny-1
            if(isinside[i, j])
                if(isinside[i-1, j] + isinside[i+1, j] + isinside[i, j-1] + isinside[i, j+1] == 0)
                    isinside[i, j] = false
                    isgrid[i, j] = false
                end
            end
        end
    end
    # Check if nodes are "chained"
    while true
        rem = Array{Float64, 2}(reshape([], 0, 2))
        for j=2:Nx-1
            for i=2:Ny-1
                if(isinside[i, j] && isinside[i+1, j] + isinside[i-1, j] + isinside[i, j+1] + isinside[i, j-1] < 2)
                    rem = vcat(rem, [i j])
                end
            end
        end
        if size(rem, 1) > 0
            for i=1:size(rem, 1)
                isinside[Int64(rem[i, 1]), Int64(rem[i, 2])] = false
                isgrid[Int64(rem[i, 1]), Int64(rem[i, 2])] = false
            end
        else
            break
        end
    end
    
    n = Array{Float64,2}(reshape([],0,4)) # defined as coordy coordx ny nx
    for j=2:Nx-1
       for i=2:Ny-1
            nx = isinside[i, j-1] - isinside[i, j+1]
            ny = isinside[i-1, j] - isinside[i+1, j]
            if((nx^2 + ny^2) > 0 && isinside[i,j] == false)
                isborder[i, j] = true
                isgrid[i, j] = true
                n = vcat(n,Array{Int64,2}([i j ny nx]))
            end
       end
    end
    
    return Nx, Ny, xrange, yrange, isinside, isborder, n, isgrid
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
    xi = sparse([], [], Float64.([]), Ny * Nx, Ny * Nx)

    #Threads.@threads for i in 1:Nx*Ny
    for k in 1:Nx*Ny
        if(isinsideflat[k])
            xi[k, k] = - (bxflat[k] +
                bxflat[k-Ny] +
                byflat[k] +
                byflat[k-1])
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
        
        xi[k, k] = - ((1.0 - nx) * bxflat[k] +
            (1.0 + nx) * bxflat[k-Ny] +
            (1.0 - ny) * byflat[k] +
            (1.0 + ny) * byflat[k-1])
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




function find_time_step(isin::BitMatrix, m::Array{Float64,2}, d::Array{Float64,2}, p::Array{Float64,2},
    bx::Array{Float64,2}, by::Array{Float64,2}, dx::Float64, alpha=0.1)
    # !!!!!!! NOT FUNCTIONAL STILL SOME WORK TO DO
    Nx = size(m,2)
    Ny = size(m,1)
    bij = zeros(size(m))
    for i = 2:Ny-2
        for j = 2:Nx-2
            bij[i,j] = bx[i,j] + bx[i,j+1] + by[i-1,j] + by[i,j]
        end
    end
    gamma = d ./ m
    println(alpha*dx^2*minimum(m[isin] ./ bij[isin]))
    println(alpha*minimum(d[isin] ./ abs.(p[isin])))
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



function get_cont_values(
    isinside::BitMatrix,
    yrange::Array{Float64, 1},
    xrange::Array{Float64, 1},
    disc_coord::Array{Float64, 2},
    disc_values::Array{Float64, 1};
    Niter = 1000
)
    
    idin = findall(isinside)
    cont_coord = zeros(length(idin), 2)
    for i = 1:length(idin)
        cont_coord[i,:] = [xrange[idin[i][2]], yrange[idin[i][1]]]
    end
    
    cont_values_vec = zeros(length(idin))
    for i = 1:size(idin, 1)
        dx = disc_coord[:,1] .- cont_coord[i,1] 
        dy = disc_coord[:,2] .- cont_coord[i,2] 
        id = argmin(dx.^2 + dy.^2)
        cont_values_vec[i] = disc_values[id]
    end
    
    cont_values = zeros(Ny, Nx)
    cont_values[idin] = cont_values_vec
    
    new_cont_values = zeros(Ny, Nx)
    tau = 0.001
    interval = 100

    @time begin
        for k in 1:Niter
            Threads.@threads for i in 2:Ny-1
                Threads.@threads for j in 2:Nx-1
                    if(isinside[i,j])
                        new_cont_values[i,j] = (1.0 - 4.0 * tau) * cont_values[i,j] + tau * (cont_values[i+1,j] +
                            cont_values[i-1,j] + cont_values[i,j+1] + cont_values[i,j-1])
                    end
                end
            end
            
            Threads.@threads for k in 1:size(n,1)
                i = Int64(n[k,1])
                j = Int64(n[k,2])
                nx = n[k,4]
                ny = n[k,3]
                if(nx == 1)
                    new_cont_values[i,j] = new_cont_values[i,j-2]
                elseif(nx == -1)
                    new_cont_values[i,j] = new_cont_values[i,j+2]
                end
                if(ny == 1)
                    new_cont_values[i,j] = new_cont_values[i-2,j]
                elseif(ny == -1)
                    new_cont_values[i,j] = new_cont_values[i+2,j]
                end
            end
            
            cont_values = copy(new_cont_values)
            
            if(mod(k,interval) == 0)
                println(k)
            end
        end
    end
    return cont_values
end
