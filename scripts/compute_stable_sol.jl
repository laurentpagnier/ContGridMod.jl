import JSON
using Plots
using HDF5
using Trapz
include("../scripts/tools.jl")
include("../scripts/get_grid.jl")
include("../scripts/get_params.jl")

border = import_border("../data/border.json")
#for dx = [40, 20, 10]
for dx = [20]
Nx, Ny, xrange, yrange, isinside, isborder, n = get_grid(border, Float64(dx))
bx, by, p, m, d = get_params(isinside, "../numerics/grid_params_" * string(dx) * ".h5")


interval = 1000
Niter = 240000
th = zeros(Ny, Nx)

@time begin
    for k in 1:Niter
        if(mod(k,interval) == 0)
            temp = copy(th)
        end
        Threads.@threads for i in 2:Ny-1
            Threads.@threads for j in 2:Nx-1
                if(isinside[i,j])
                    bij = (by[i-1,j] + by[i,j] + bx[i,j] + bx[i,j+1])
                    th[i,j] = (by[i,j] * th[i+1,j] + by[i-1,j] * th[i-1,j] + 
                        bx[i,j+1] * th[i,j+1] + bx[i,j] * th[i,j-1] + dx^2*p[i,j]) / bij
                end
            end
        end
        Threads.@threads for k in 1:size(n,1)
        #for k in 1:size(n,1)
            i = Int64(n[k,1])
            j = Int64(n[k,2])
            nx = n[k,4]
            ny = n[k,3]
            if(nx == 1)
                th[i,j] = th[i,j-2]
            elseif(nx == -1)
                th[i,j] = th[i,j+2]
            end
            if(ny == 1)
                th[i,j] = th[i-2,j]
            elseif(ny == -1)
                th[i,j] = th[i+2,j]
            end
        end  
        if(mod(k,interval) == 0)
            println([k maximum(abs.(th-temp))])
        end

    end
end

fid = h5open("../numerics/stable_" * string(dx) * ".h5", "w")
write(fid, "th", th)
close(fid)
end

