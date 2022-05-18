all:
	@echo select an option:
	@echo '"make grid" to compute grid parameters'
	@echo '"make stable" to compute stable solutions'
install:
	@test -e numerics || mkdir numerics
clear:
	rm numerics/*
grid: num
	cd scripts;\
	julia -t $$(nproc) compute_grid_params.jl
stable:
	cd scripts;\
	julia -t $$(nproc) compute_stable_sol.jl
dyn:
	cd scripts;\
	julia -t $$(nproc) compute_dynamics.jl
download:
	wget 
