export add_local_disturbance!, add_local_disturbance_with_gps!


function add_local_disturbance!(
    contmod::ContModel,
    coord::Tuple{Float64, Float64},
    dP::Float64;
    tstart::Float64 = 0.0,
)
    id = find_nearest(coord, contmod.mesh.node_coord)
    temp = zeros(contmod.mesh.Nnode)
    temp[id] = dP
    dp(t) = temp .* (t .> tstart)
    contmod.dp = dp
    nothing
end


function add_local_disturbance_with_gps!(
    contmod::ContModel,
    gps_coord::Tuple{Float64, Float64},
    dP::Float64,
)
    c = albers_projection([gps_coord[1]  gps_coord[2]] ./ (180 / pi) )
    c ./= contmod.mesh.scale_factor
    add_local_disturbance!(contmod, (c[1], c[2]), dP)
    nothing
end
