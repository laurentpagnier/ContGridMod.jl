function local_disturbance(
    grid_coord::Array{Float64, 2},
    location::Array{Float64,1},
    dP::Float64,
    sigma::Float64,
    dx::Float64
)
    dp = zeros(size(grid_coord,1))
    for k = 1:size(grid_coord,1)
        dp[k] = exp(-((grid_coord[k,1] - location[1])^2 +(grid_coord[k,2] - location[2])^2)/2/sigma^2)
    end
    return dP .* dp ./ sum(dp) / dx^2
end
