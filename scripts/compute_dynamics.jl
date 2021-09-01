# load libraries and scripts
import JSON
using Plots
using HDF5
using Trapz
include("../scripts/tools.jl")
include("../scripts/get_grid.jl")
include("../scripts/get_params.jl")
include("../scripts/disturbances.jl")

#load borders
border = import_border("../data/border.json")

# create the lattice grid
dx = 20
Nx, Ny, xrange, yrange, isinside, isborder, n = get_grid(border, Float64(dx))
bx, by, p, m, d = get_params(isinside, "../numerics/grid_params_" * string(dx) * ".h5")

# define a disturbance
#dP = -9.0
dP = 0.0
sigma = 30.0
location = [-1500., -900.]
dp = local_disturbance(isinside, xrange, yrange, location, dP, sigma)
do_plot(isinside, dp)

# perform a dynamical simulation
interval = 200
dt = 0.0001
Ndt = 20000

println(maximum(u)*dt^2/dx^2)

#gamma = 100*d ./ m
#m = m ./10

omegas = zeros(Ny,Nx,1 + Int64(ceil(Ndt/interval)))
thetas = zeros(Ny,Nx,1 + Int64(ceil(Ndt/interval)))
th_new = zeros(Ny,Nx)
reset = true
data = h5read("../numerics/stable_" * string(dx) * ".h5", "/")
if(reset)
    th_old = copy(data["th"])
    th = copy(data["th"])
    omegas[:,:,1] = zeros(size(th))
    thetas[:,:,1] = copy(th)  
else
    omegas[:,:,1] = (th - th_old) / dt
    thetas[:,:,1] = copy(th)   
end

ts = zeros(1 + Int64(ceil(Ndt/interval)))
chi = 1 .+ gamma*dt/2

@time begin
    for t in 1:Ndt
        Threads.@threads for i in 2:Ny-1
            Threads.@threads for j in 2:Nx-1
                if(isinside[i,j])
                    bij = (by[i-1,j] + by[i,j] + bx[i,j] + bx[i,j+1])
                    th_new[i,j] = (2 - bij / m[i,j] * dt^2 / dx^2) / chi[i,j] * th[i,j] + 
                        (gamma[i,j] * dt/2 - 1) /chi[i,j] * th_old[i,j] + 
                        dt^2 / dx^2 / chi[i,j] / m[i,j] * 
                        (by[i,j] * th[i+1,j] + by[i-1,j] * th[i-1,j] +
                        bx[i,j+1] * th[i,j+1] + bx[i,j] * th[i,j-1]) +
                        dt^2 / chi[i,j] / m[i,j] * (p[i,j] + dp[i,j])
                end
            end
        end   

        # impose boundary condition
        Threads.@threads for k in 1:size(n,1)
            i = Int64(n[k,1])
            j = Int64(n[k,2])
            nx = n[k,4]
            ny = n[k,3]
            if(nx == 1)
                th_new[i,j] = th_new[i,j-2]
            else(nx == -1)
                th_new[i,j] = th_new[i,j+2]
            end
            if(ny == 1)
                th_new[i,j] = th_new[i-2,j]
            else(ny == -1)
                th_new[i,j] = th_new[i+2,j]
            end
        end
        
        if(mod(t,interval) == 0)
            println("NIter: ", t)
            omegas[:,:,Int64(t/interval) + 1] = (th_new-th) / dt
            thetas[:,:,Int64(t/interval) + 1] = th_new
            ts[Int64(t/interval) + 1] = t*dt
        end
        th_old = copy(th)
        th = copy(th_new)
    end
end


