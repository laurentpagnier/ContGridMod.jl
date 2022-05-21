export add_local_disturbance!, add_local_disturbance_with_gps!


function add_local_disturbance!(
    contmod::ContModel,
    coord::Array{Float64,1},
    dP::Float64,
    sigma::Float64,
)
    grid_coord = contmod.mesh.coord
    dp = zeros(size(grid_coord,1))
    for k = 1:size(grid_coord,1)
        dp[k] = exp(-((grid_coord[k,1] - coord[1])^2 + (grid_coord[k,2] - coord[2])^2) / 2 / sigma^2)
    end
    
    contmod.p += dP .* dp ./ sum(dp) / contmod.mesh.dx^2
end


function add_local_disturbance_with_gps!(
    contmod::ContModel,
    gps_coord::Array{Float64,1},
    dP::Float64,
    sigma::Float64,
)
    grid_coord = contmod.mesh.coord
    coord = albers_projection(reshape(gps_coord,1,2) ./ (180 / pi) )
    coord ./= contmod.scale_factor
    dp = zeros(size(grid_coord,1))
    for k = 1:size(grid_coord,1)
        dp[k] = exp(-((grid_coord[k,1] - coord[1])^2 + (grid_coord[k,2] - coord[2])^2) / 2 /sigma^2)
    end
    
    contmod.p += dP .* dp ./ sum(dp) / contmod.mesh.dx^2
end
