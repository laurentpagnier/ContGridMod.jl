export get_mesh

function get_mesh(
    border::Array{Float64, 2},
    dx::Float64
)
    xlim = [minimum(border[:, 1])-2*dx maximum(border[:, 1])+2*dx]
    ylim = [minimum(border[:, 2])-2*dx maximum(border[:, 2])+2*dx]
    xrange = collect(xlim[1]:dx:xlim[2])
    yrange = collect(ylim[1]:dx:ylim[2])

    Nx = length(xrange)
    Ny = length(yrange)

    isborder = falses(Ny, Nx)
    isinside = falses(Ny, Nx)
    isgrid = falses(Ny, Nx)

    # check if inside the shape
    for j = 1:Nx
        for i = 1:Ny
            p = [xrange[j] yrange[i]]
            isinside[i, j] = inPolygon(p, border)[1]
        end
    end
    
    # check if alone 
    for j=2:Nx-1
        for i=2:Ny-1
            if(isinside[i, j])
                if(isinside[i-1, j] + isinside[i+1, j] + isinside[i, j-1] + isinside[i, j+1] == 0)
                    isinside[i, j] = false
                end
            end
        end
    end
    
    # add a boundary layer encompassing the inside
    for j=2:Nx-1
       for i=2:Ny-1
            nx = isinside[i, j-1] - isinside[i, j+1]
            ny = isinside[i-1, j] - isinside[i+1, j]
            if((nx^2 + ny^2) > 0 && isinside[i,j] == false)
                isgrid[i, j] = true
            end
       end
    end
    
    isgrid = isgrid .| isinside
    # check if inside and considered as outside. This check isnt able
    # to find hole that are bigger than a single point. Som it should
    # be improved in the future
    for j=2:Nx-1
        for i=2:Ny-1
            if(.!isgrid[i, j])
                if(isgrid[i-1, j] + isgrid[i+1, j] + isgrid[i, j-1] + isgrid[i, j+1] == 4)
                    isgrid[i, j] = true
                    isinside[i, j] = true
                end
            end
        end
    end
    
    for j=2:Nx-1
        for i=2:Ny-1
            if(isgrid[i, j])
                if(!isgrid[i-1, j] + !isgrid[i+1, j] + !isgrid[i, j-1] + !isgrid[i, j+1] > 2)
                    isgrid[i, j] = false
                end
            end
        end
    end
    
    n = Array{Float64, 2}( reshape([], 0, 4) ) # defined as coordy coordx ny nx
    for j=2:Nx-1
        for i=2:Ny-1
            if(isgrid[i, j])
                if(!isgrid[i,j+1] & isgrid[i,j-1])
                    nx = 1
                elseif(!isgrid[i,j-1] & isgrid[i,j+1])
                    nx = -1
                else
                    nx = 0
                end
                if(!isgrid[i+1,j] & isgrid[i-1,j])
                    ny = 1
                elseif(!isgrid[i-1,j] & isgrid[i+1,j])
                    ny = -1
                else
                    ny = 0
                end
                if((nx^2 + ny^2) > 0)
                    isborder[i, j] = true
                    n = vcat(n, Array{Int64,2}([i j ny nx]))
                end
            end
        end
    end

    coord = reshape([],0,2)
    for j in 1:Nx
        for i in 1:Ny
            coord = [coord; reshape([yrange[i] xrange[j]], 1, 2)]
        end
    end

    isinside = isgrid .& .!isborder
    
    mesh = Mesh(Nx, Ny, Float64.(coord), vec(isinside),
        vec(isborder), vec(isgrid), xrange, yrange, n, dx)
    return mesh
end
