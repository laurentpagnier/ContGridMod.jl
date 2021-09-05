function set_ref_phase(
    isinside::BitMatrix,
    yrange::Array{Float64, 1},
    xrange::Array{Float64, 1},
    theta::Array{Float64, 2},
    coord::Array{Float64, 1},
    th_ref::Float64 = 0.0
)
    # set the phases according to some reference phase at at a given location 
    idin = findall(isinside)
    cont_coord = zeros(length(idin), 2)
    for i = 1:length(idin)
        cont_coord[i, :] = [xrange[idin[i][2]], yrange[idin[i][1]]]
    end
    dx = cont_coord[:, 1] .- coord[1]
    dy = cont_coord[:, 2] .- coord[2]
    id = argmin(dx.^2 + dy.^2)
    println(cont_coord[id,:], theta[idin[id]])
    theta[isinside] .= theta[isinside] .- theta[idin[id]] .+ th_ref
    return theta
end


function power_flow(
    isinside::BitMatrix,
    theta::Array{Float64, 2},
    bx::Array{Float64, 2},
    by::Array{Float64, 2},
    dx::Float64
)
    # compute the power flow vector field in the continuous model
    dx_th = zeros(size(theta))
    dy_th = zeros(size(theta))
    for i in 2:size(theta, 1)-1
        for j in 2:size(theta, 2)-1
            if(isinside[i, j])
                dx_th[i, j] = (theta[i, j+1] - theta[i ,j-1]) / 2.0 / dx
                dx_th[i, j] = (theta[i, j+1] - theta[i ,j-1]) / 2.0 / dx
                dy_th[i, j] = (theta[i+1, j] - theta[i-1, j]) / 2.0 / dx
            end
        end
    end
    fx = bx .* dx_th
    fy = by .* dy_th
    return fx, fy
end


function get_discrete_values(
    isinside::BitMatrix,
    yrange::Array{Float64, 1},
    xrange::Array{Float64, 1},
    cont_values::Array{Float64, 2},
    disc_coord::Array{Float64, 2}
)
    
    idin = findall(isinside)
    cont_coord = zeros(length(idin), 2)
    for i = 1:length(idin)
        cont_coord[i,:] = [xrange[idin[i][2]], yrange[idin[i][1]]]
    end
    
    disc_values = zeros(size(coord, 1))
    for i = 1:size(disc_coord, 1)
        dx = cont_coord[:, 1] .- disc_coord[i, 1]
        dy = cont_coord[:, 2] .- disc_coord[i, 2]
        dist = min.(
            sqrt.(dx.^2 + dy.^2),
            10*minimum(sqrt.(dx.^2 + dy.^2))
            ) # max is here to prevent NaN
        factor = exp.(-dist) ./ sum(exp.(-dist)) # here we use softmax for no particular reason
        if(isnan(sum(factor)))
            println(minimum(dist))
            println(maximum(dist))
        end
        disc_values[i] = sum(factor .* cont_values[idin])
    end
    
    return disc_values
end
