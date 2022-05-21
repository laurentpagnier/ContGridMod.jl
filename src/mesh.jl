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
    
    # create incidence matrix
    bus_id = Int64.(zeros(0,2))
    bus_coord = zeros(0,2)
    for j=1:Nx-1
        for i=1:Ny-1
            if(isgrid[i, j])
                bus_id = [bus_id; j i]
                bus_coord = [bus_coord; yrange[i] xrange[j]]
            end
        end
    end

    incidence_mat = Int64.(zeros(0,2))
    line_coord = zeros(0,4)
    for j=1:Nx-1
        for i=1:Ny-1
            if(isgrid[i, j] && isgrid[i, j+1])
                id1 = findfirst(all(bus_id .== [j i],dims=2))[1]
                id2 = findfirst(all(bus_id .== [j+1 i],dims=2))[1]
                incidence_mat = [incidence_mat; id1 id2]
                line_coord = [line_coord; [yrange[i] xrange[j] yrange[i] xrange[j+1]]]
            end
            if(isgrid[i, j] && isgrid[i+1, j])
                id1 = findfirst(all(bus_id .== [j i],dims=2))[1]
                id2 = findfirst(all(bus_id .== [j i+1],dims=2))[1]
                incidence_mat = [incidence_mat; id1 id2]
                line_coord = [line_coord; [yrange[i] xrange[j] yrange[i+1] xrange[j]]]
            end
        end
    end
    
    return mesh = Mesh(
        Nx,
        Ny,
        Float64.(bus_coord[:,2:-1:1]),
        Float64.(line_coord),
        incidence_mat,
        vec(isgrid),
        yrange,
        xrange,
        dx,
    )
end
