import JSON
using Plots
using HDF5
using Trapz
include("tools.jl")
include("get_grid.jl")
include("get_params.jl")

sigma = 100.

#for dx = [40, 20, 10]
for dx = [20]
	border = import_border("../data/border.json")
	Nx, Ny, xrange, yrange, isinside, isborder, n = get_grid(border, Float64(dx))
	bx, by, p, m, d = get_params(isinside, isborder, sigma, Float64(dx), yrange, xrange, "../data/pantagruel.h5",
        	"../numerics/grid_params_" * string(dx) * ".h5")
end
