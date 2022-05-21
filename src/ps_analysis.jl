export set_ref_phase!

function set_ref_phase!(
    contmod::ContModel,
    coord_ref::Array{Float64, 1};
    th_ref::Float64 = 0.0
)
    # set the phases according to some reference phase at at a given location
    id = argmin((contmod.mesh.coord[:, 1] .- coord_ref[1]).^2 +
        (contmod.mesh.coord[:, 2] .- coord_ref[2]).^2)
    n = length(contmod.th)
    M = sparse([1:n; 1:n], [1:n; id*ones(n)], [ones(n); -ones(n)])
    contmod.th = M * contmod.th + th_ref * ones(n)
end


function set_ref_phase_with_gps!(
    contmod::ContModel,
    gps_coord_ref::Array{Float64, 1};
    th_ref::Float64 = 0.0
)
    # set the phases according to some reference phase at at a given location
    gps_coord_ref = reshape(gps_coord_ref, 1, 2)
    id = find_node_from_gps(contmod, gps_coord_ref)[1]
    n = length(contmod.th)
    M = sparse([1:n; 1:n], [1:n; id*ones(n)], [ones(n); -ones(n)])
    contmod.th = M * contmod.th + th_ref * ones(n)
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
    
    temp = findall(isinside)
    idin = Int64.(zeros(size(temp, 1), 2))
    for i in 1:size(temp,1)
        idin[i,:] = [temp[i][1], temp[i][2]]
    end
    
    Threads.@threads for k in 1:size(idin, 1)
        i = idin[k, 1]
        j = idin[k, 2]
        if(abs(theta[i,j]) < 1E-6)
            println([i,j])
        end
        #dx_th[i, j] = (theta[i, j] - theta[i, j-1]) / dx
        #dy_th[i, j] = (theta[i, j] - theta[i-1, j]) / dx
        dx_th[i, j] = (theta[i, j+1] - theta[i ,j-1]) / 2.0 / dx
        dy_th[i, j] = (theta[i+1, j] - theta[i-1, j]) / 2.0 / dx
    end
    
    dx_th_avg = sum(dx_th[isinside]) / sum(isinside)
    #dx_th[dx_th .< 1 * dx_th_avg] .= 0
    #dy_th[dy_th .< 1 * dx_th_avg] .= 0
    
    #for i in 2:size(theta, 1)-1
    #    for j in 2:size(theta, 2)-1
    #        if(isinside[i, j])
    #            #dx_th[i, j] = (theta[i, j+1] - theta[i ,j-1]) / 2.0 / dx
    #            #dy_th[i, j] = (theta[i+1, j] - theta[i-1, j]) / 2.0 / dx
    #            dx_th[i, j] = (theta[i, j] - theta[i ,j-1]) / dx
    #            dy_th[i, j] = (theta[i, j] - theta[i-1, j]) / dx
    #        end
    #    end
    #end
    fx = bx .* dx_th
    fy = by .* dy_th
    
    #fx = dx_th
    #fy = dy_th
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
