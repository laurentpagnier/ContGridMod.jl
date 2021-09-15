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

