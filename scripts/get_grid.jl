function get_grid(border::Array{Float64,2}, dx::Float64)

    xlim = [minimum(border[:,1])-dx maximum(border[:,1])+dx]
    ylim = [minimum(border[:,2])-dx maximum(border[:,2])+dx]
    xrange = collect(xlim[1]:dx:xlim[2])
    yrange = collect(ylim[1]:dx:ylim[2])

    Nx = length(xrange)
    Ny = length(yrange)

    isgrid = Bool.(zeros(Ny,Nx))
    isborder = Bool.(zeros(Ny,Nx))

    for j=1:Nx
        for i=1:Ny
            p = [xrange[j] yrange[i]]
            isgrid[i,j] = inPolygon(p,border)[1]
        end
    end

    # boundary normal vector
    n = Array{Float64,2}(reshape([],0,4)) # defined as coordy coordx nx ny
    for j=2:Nx-1
       for i=2:Ny-1
            if((sum(isgrid[i+1,j] + isgrid[i-1,j] + isgrid[i,j+1] + isgrid[i,j-1]) < 4) & isgrid[i,j])
                isborder[i,j] = 1
                nx = isgrid[i,j-1]-isgrid[i,j+1]
                ny = isgrid[i-1,j]-isgrid[i+1,j]
                n = vcat(n,Array{Float64,2}([i j nx ny]))
            end
       end
    end
    return Nx, Ny, xrange, yrange, isgrid, isborder, n
end

