# load libraries and scripts
import JSON
using Plots
using HDF5
using Trapz
using Statistics
using SparseArrays
using LinearAlgebra
include("../scripts/tools.jl")
include("../scripts/get_grid.jl")
include("../scripts/get_params.jl")
include("../scripts/disturbances.jl")
include("../scripts/stable.jl")
include("../scripts/dynamics.jl")
include("../scripts/vectorize.jl")

#load borders
border = import_border("../data/borders/border.json");

# create the lattice grid
dx = 20.0
Nx, Ny, xrange, yrange, isinside, isborder, n, isgrid = get_grid(border, Float64(dx));


m = 1e-5 * ones(Ny, Nx)
d = 0.3 * m
bx = 8 * ones(Ny, Nx)
by = 8 * ones(Ny, Nx)
p = zeros(Ny, Nx)
m[.!isgrid] .= 0
d[.!isgrid] .= 0
p[.!isgrid] .= 0;

isinsideflat, pflat, minvflat, gammaflat, xi = vectorize_ref(isinside, isborder, n, bx, by, p, m, d);
isgridflat = vec(isinside .| isborder);

# define a disturbance
dP = -9.0
# dP = 0.0
sigma = 100.0
location = [-1500.0, -900.0]
dp = local_disturbance(isgrid, xrange, yrange, location, dP, sigma)
dpflat = vec(dp)
println("Synchronized frequency: ", trapz((yrange, xrange), p .+ dp) / trapz((yrange, xrange), d))

th0 = zeros(Ny*Nx)
#ts_be, ~, omegas_be = perform_dyn_sim_vec_backward_euler(isgridflat, xi, pflat[isgridflat]+dpflat[isgridflat],
#    minvflat, gammaflat, th0, interval = 5, Ndt = 2500, dt = 0.01)
ts_cn, ~, omegas_cn = perform_dyn_sim_vec_crank_nicolson(isgridflat, xi, pflat[isgridflat]+dpflat[isgridflat],
    minvflat, gammaflat, th0, interval = 5, Ndt = 500, dt = 0.05)
#ts_v, ~, omegas_v = perform_dyn_sim_vec(isgridflat, xi, pflat[isgridflat]+dpflat[isgridflat],
#    minvflat, gammaflat, th0, interval = 1000, Ndt = 250000, dt = 0.0001)
#omegas_be = back_to_2d(isgrid, omegas_be)




