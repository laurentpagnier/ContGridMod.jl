function get_grid(
    border::Array{Float64, 2},
    dx::Float64
)
    xlim = [minimum(border[:, 1])-2*dx maximum(border[:, 1])+2*dx]
    ylim = [minimum(border[:, 2])-2*dx maximum(border[:, 2])+2*dx]
    xrange = collect(xlim[1]:dx:xlim[2])
    yrange = collect(ylim[1]:dx:ylim[2])

    Nx = length(xrange)
    Ny = length(yrange)

    isborder = Bool.(zeros(Ny,Nx))
    isinside = Bool.(zeros(Ny,Nx))

    for j = 1:Nx
        for i = 1:Ny
            p = [xrange[j] yrange[i]]
            isinside[i, j] = inPolygon(p, border)[1]
        end
    end

    n = Array{Float64, 2}( reshape([], 0, 4) ) # defined as coordy coordx ny nx
    for j=2:Nx-1
        for i=2:Ny-1
            if(isinside[i,j] & !isinside[i,j+1] & isinside[i,j-1])
                nx = 1
            elseif(isinside[i,j] & !isinside[i,j-1] & isinside[i,j+1])
                nx = -1
            else
                nx = 0
            end
            if(isinside[i,j] & !isinside[i+1,j] & isinside[i-1,j])
                ny = 1
            elseif(isinside[i,j] & !isinside[i-1,j] & isinside[i+1,j])
                ny = -1
            else
                ny = 0
            end
            if((nx^2 + ny^2) > 0)
                if((isinside[i-1,j] + isinside[i+1,j] + isinside[i,j-1] + isinside[i,j+1]) < 2)
                    isinside[i,j] = false
                    if(nx^2 > 0)
                        isborder[i, j-Int64(nx)] = true
                        n = vcat(n, Array{Int64,2}([i j-Int64(nx) ny nx]))
                    end
                    if(ny^2 > 0)
                        isborder[i-Int64(ny), j] = true
                        n = vcat(n, Array{Int64,2}([i-Int64(ny) j ny nx]))
                    end
                else
                    isborder[i, j] = true
                    n = vcat(n, Array{Int64,2}([i j ny nx]))
                end
            end
        end
    end
    
    isinside = isinside .& .!isborder
    return Nx, Ny, xrange, yrange, isinside, isborder, n
end

