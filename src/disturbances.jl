function add_local_disturbance!(
    contmod::ContModel,
    location::Array{Float64,1},
    dP::Float64,
    sigma::Float64,
)
    grid_coord = contmod.coord[contmod.isgrid,:]
    dp = zeros(size(grid_coord,1))
    for k = 1:size(grid_coord,1)
        dp[k] = exp(-((grid_coord[k,1] - location[1])^2 + (grid_coord[k,2] - location[2])^2)/2/sigma^2)
    end
    
    contmod.p += dP .* dp ./ sum(dp) / contmod.dx^2
end
