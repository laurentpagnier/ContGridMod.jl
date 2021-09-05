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
