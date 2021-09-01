all:
	@echo select an option:
	@echo '"make grid" to compute grid parameters'
	@echo '"make stable" to compute stable solutions'
num:
	@test -e numerics || mkdir numerics
clear:
	rm numerics/*
grid: num
	cd scripts;\
	julia -t 8 compute_grid_params.jl
stable:
	cd scripts;\
	julia -t 8 compute_stable_sol.jl
dyn:
	cd scripts;\
	julia -t 8 compute_dynamics.jl
