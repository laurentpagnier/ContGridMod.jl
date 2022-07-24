export back_to_2d, albers_projection, import_json_numerics, import_border, get_discrete_values, copy_model

function inPolygon(
    p::Matrix{Float64},
    poly::Matrix{Float64},
)
    N = size(poly, 1)
    isin = falses(size(p, 1))
    for k = 1:size(p, 1)
        j = N
        for i = 1:N
            if (
                ((poly[i, 2] < p[k, 2]) & (poly[j, 2] >= p[k, 2])) |
                ((poly[j, 2] < p[k, 2]) & (poly[i, 2] >= p[k, 2]))
            )
                if (
                    poly[i, 1] +
                    (p[k, 2] - poly[i, 2]) / (poly[j, 2] - poly[i, 2]) *
                    (poly[j, 1] - poly[i, 1]) < p[k, 1]
                )
                    isin[k] = !isin[k]
                end
            end
            j = i
        end
    end
    return isin
end


function copy_model(
    cm::ContModel
)
    keys = fieldnames(ContGridMod.ContModel)
    data = Tuple(perform_copy(getfield(cm, k)) for k in keys)
    return ContModel(data...)
end


function copy_mesh(
    mesh::Mesh
)
    keys = fieldnames(Mesh)
    data = Tuple(copy(getfield(mesh, k)) for k in keys)
    return Mesh(data...)
end


function perform_copy(field)
    if typeof(field) == Mesh
        return copy_mesh(field)
    else
        return copy(field)
    end
end


function albers_projection(
    coord::Array{Float64,2}; # as latitude, longitude 
    lon0::Float64 = 13.37616 / 180 * pi,
    lat0::Float64 = 46.94653 / 180 * pi,
    lat1::Float64 = 10 / 180 * pi,
    lat2::Float64 = 50 / 180 * pi,
    R::Float64 = 6371.0
)
    # see https://en.wikipedia.org/wiki/Albers_projection
    n = 1 / 2 * (sin(lat1) + sin(lat2))
    theta = n * (coord[:, 2] .- lon0)
    c = cos(lat1)^2 + 2 * n * sin(lat1)
    rho = R / n * sqrt.(c .- 2 * n * sin.(coord[:, 1]))
    rho0 = R / n * sqrt.(c - 2 * n * sin(lat0))
    x = rho .* sin.(theta)
    y = rho0 .- rho .* cos.(theta)
    return [x y]
end


function import_border(
    filename::String
)
    data = JSON.parsefile(filename)
    N = size(data["border"], 1)

    b = zeros(N, 2)
    for i = 1:N
        b[i, 1] = data["border"][i][1] / 180 * pi
        b[i, 2] = data["border"][i][2] / 180 * pi
    end

    if (b[1, :] != b[end, :])
        b = vcat(b, reshape(b[1, :], 1, 2))
    end

    b = albers_projection(b)
    scale_factor = max(
        maximum(b[:,1]) - minimum(b[:,1]),
        maximum(b[:,2]) - minimum(b[:,2])
        )
    return b / scale_factor, scale_factor
end


function get_discrete_id(
    grid_coord::Array{Float64,2},
    disc_coord::Array{Float64,2},
)
    ids = Int64.(zeros(size(disc_coord,1)))
    for i in 1:size(ids,1)
        ids[i] = argmin((grid_coord[:,1] .- disc_coord[i,1]).^2 +
        (grid_coord[:,2] .- disc_coord[i,2]).^2)
    end

    return ids
end


function get_discrete_values(
    ids::Vector{Int64},
    v::Vector{Float64}
)
    return [v[ids[i]] for i = 1:length(ids)]
end


function get_discrete_values(
    grid_coord::Matrix{Float64},
    disc_coord::Matrix{Float64},
    v::Vector{Float64}
)
    disc_v = zeros(size(disc_coord,1))
    #map((i) -> )
    for i in 1:size(disc_v,1)
        id = argmin((grid_coord[:,1] .- disc_coord[i,1]).^2 +
        (grid_coord[:,2] .- disc_coord[i,2]).^2)
        disc_v[i] = disc_v[i] + v[id] # changed to avoid "Mutating arrays is not supported"
    end
    #=
    for i in 1:size(disc_v,1)
        id = argmin((grid_coord[:,1] .- disc_coord[i,1]).^2 +
        (grid_coord[:,2] .- disc_coord[i,2]).^2)
        disc_v[i] = disc_v[i] + v[id] # changed to avoid "Mutating arrays is not supported"
    end
    =#
    return disc_v
end


function import_json_numerics(
    filename::String
)
    raw_data = JSON.parsefile(filename)
    keys = raw_data.keys
    isdef = falses(length(keys))
    for i = 1:length(keys)
        isdef[i] = isassigned(keys, i)
    end
    keys = keys[isdef]
    data = Dict{String, Array{Float64,2}}()
    for k  in keys
        N = size(raw_data[k], 1)
        M = size(raw_data[k][1], 1) # assumes that data consists of matrices
        temp = zeros(N, M)
        for i = 1:N
            for j = 1:M
                temp[i, j] = raw_data[k][i][j]
            end
        end
        data[k] = temp
    end
    return data, keys
end


function back_to_2d(
    isgrid::BitVector,
    Ny::Int64,
    Nx::Int64,
    valueflat::Matrix{Float64},
)
    value = zeros(Ny, Nx, size(valueflat, 2))
    for t in 1:size(valueflat, 2)
    	temp = zeros(size(isgrid))
    	temp[isgrid] .= valueflat[:, t]
        value[:,:, t] = reshape(temp, Ny, Nx)
    end
    return value
end


function get_cont_values(
    isgrid::BitVector,
    grid_coord::Array{Float64, 2},
    disc_coord::Array{Float64, 2},
    disc_values::Array{Float64, 1}
)
    v = zeros(size(grid_coord,1))
    for i = 1:size(grid_coord, 1)
        id = argmin((grid_coord[i,1] .- disc_coord[:,1]).^2 +
            (grid_coord[i,2] .- disc_coord[:,2]).^2)
        v[i] = disc_values[id]
    end
    values = zeros(size(isgrid,1))
    values[isgrid] = v
    return values
end


function find_gen(
    dm::DiscModel,
    gps_coord::Array{Float64, 2},
    dP::Float64;
    scale_factor::Float64 = 1.0
)
    #find the the nearest generator that can "withstand" a dP fault
    coord = albers_projection(gps_coord[:,[2;1]] ./ (180 / pi) )
    coord = coord[:,[2,1]] / scale_factor
    idprod = findall((dm.gen .> 0.0))
    idavail = findall(dm.max_gen[idprod] .> dP) # find id large gens in the prod list
    # println(size(idprod))
    # println(idavail)
    #idprod = findall((dm.gen .> 0.0))
    id = Int64.(zeros(size(coord,1)))  # index in the full gen list
    id2 = Int64.(zeros(size(coord,1))) # index in the producing gen list
    for i in 1:size(coord,1)
        temp = dm.idgen[idprod[idavail]]
        id[i] = idprod[idavail[argmin((dm.coord[temp,1] .- coord[i,1]).^2 +
            (dm.coord[temp,2] .- coord[i,2]).^2)]]
        id2[i] = idavail[argmin((dm.coord[temp,1] .- coord[i,1]).^2 +
            (dm.coord[temp,2] .- coord[i,2]).^2)]
    end
    return id, id2
end


function find_node_from_gps(
    cm::ContModel,
    gps_coord::Array{Float64, 2}
)
    coord = albers_projection(gps_coord ./ (180 / pi) )
    coord ./= cm.mesh.scale_factor

    id = Int64.(zeros(size(coord,1))) # index in the full gen list
    for i in 1:size(coord,1)
        id[i] = argmin((cm.mesh.coord[:,1] .- coord[i,1]).^2 +
            (cm.mesh.coord[:,2] .- coord[i,2]).^2)
    end
    return id
end


function find_node(
    cm::ContModel,
    coord::Array{Float64, 2}
)
    id = Int64.(zeros(size(coord,1))) # index in the full gen list
    for i in 1:size(coord,1)
        id[i] = argmin((cm.mesh.coord[:,1] .- coord[i,1]).^2 +
            (cm.mesh.coord[:,2] .- coord[i,2]).^2)
    end
    return id
end


function find_node(
    m::Mesh,
    coord::Array{Float64, 2}
)
    id = Int64.(zeros(size(coord,1))) # index in the full gen list
    for i in 1:size(coord,1)
        id[i] = argmin((m.coord[:,1] .- coord[i,1]).^2 +
            (m.coord[:,2] .- coord[i,2]).^2)
    end
    return id
end

