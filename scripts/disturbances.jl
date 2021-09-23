function local_disturbance_fast(
    grid_coord::Array{Tuple{Float64, Float64},1},
    location::Array{Float64,1},
    dP::Float64,
    sigma::Float64,
    dx::Float64
)
    dp = zeros(length(grid_coord))
    for k = 1:length(grid)
        dp[k] = exp(-((grid_coord[i][1] - location[1])^2 +(grid_coord[i][2] - location[2])^2)/2/sigma^2)
    end
    return dP .* dp ./ sum(dp) * dx^2
end



function local_disturbance(
    isinside::BitMatrix,
    xrange::Array{Float64,1},
    yrange::Array{Float64,1},
    location::Array{Float64,1},
    dP::Float64,
    sigma::Float64
)
    dp = zeros(Ny, Nx)
    for i = 1:Nx
        for j = 1:Ny
            if(isinside[j,i])
                dp[j, i] = exp(-((xrange[i] - location[1])^2 +(yrange[j] - location[2])^2)/2/sigma^2)
            end
        end
    end
    return dP .* dp ./ trapz((yrange, xrange), dp)
end
