export back_to_2d, albers_projection, import_json_numerics, import_border,
    get_discrete_values, copy_model, to_mat, load_discrete_model

using HDF5

function in_polygon(
    point::Vector{Tuple{Float64, Float64}},
    border_point::Vector{Tuple{Float64, Float64}},
)
    is_in = Bool[]
    for p in point
        temp = false
        prev = border_point[end]
        for b = border_point
            if (
                ((b[2] < p[2]) & (prev[2] >= p[2])) |
                ((prev[2] < p[2]) & (b[2] >= p[2]))
            )
                if (
                    b[1] +
                    (p[2] - b[2]) / (prev[2] - b[2]) *
                    (prev[1] - b[1]) < p[1]
                )
                    temp = !temp
                end
            end
            prev = b
        end
        push!(is_in, temp)
    end
    return is_in
end

to_mat(v) = mapreduce(x->[x[1], x[2]]',vcat, v)


function select_observation_points(dm, mesh)
    id_obs = Int64[]
    id_obs_cm = Int64[]
    for (i, v) in enumerate(mesh.cell_vertices)
        list = findall(ContGridMod.in_polygon(dm.coord, v))
        if !isempty(list)
            push!(id_obs, rand(list))
            push!(id_obs_cm, i)
        end
    end
    return id_obs, id_obs_cm
end


function find_nearest(
    p::Tuple{Float64,Float64},
    ref::Vector{Tuple{Float64,Float64}},
)
    d_min, id = Inf, 0
    for (i, r) in enumerate(ref)
        d = (r[1] - p[1])^2 + (r[2] - p[2])^2
        if d < d_min
            d_min, id = d, i
        end
    end
    return id
end


function find_node_from_gps(
    contmod::ContModel,
    gps_coord::Vector{Tuple{Float64, Float64}},
)
    coord = albers_projection(to_mat(gps_coord) ./ (180 / pi) )
    coord = eachrow(coord) .|> c -> (c[1], c[2]) ./ contmod.mesh.scale_factor
    return coord .|> c-> find_nearest(c, contmod.mesh.node_coord)
end


function albers_projection(
    coord::Matrix{Float64}; # as latitude, longitude 
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
    filename::String,
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


function load_discrete_model(
    dataname::String,
    scaling_factor::Float64
)
    d = h5read(dataname, "/")
    line_list = eachrow(d["line_list"]) .|> r -> (r[1], r[2]) 
    coord = ContGridMod.albers_projection(d["coord"] ./ (180 / pi))
    coord = eachrow(coord) .|> r -> (r[1], r[2]) ./ scaling_factor

    dm = ContGridMod.DiscModel(d["m_gen"], d["d_gen"], d["id_gen"], d["id_slack"], coord,
        d["d_load"], line_list, d["b"], vec(d["p_load"]), d["th"],
        d["p_gen"], d["max_gen"], d["Nbus"], d["Ngen"], d["Nline"])
    return dm
end


function load_matlab_model(
    dataname::String,
    scaling_factor::Float64
)
    data = h5read(dataname, "/")
    coord = albers_projection( data["bus_coord"] ./ (180 / pi) )
    coord = eachrow(coord) .|> c -> (c[1], c[2]) ./ scaling_factor
    line_list = eachrow(data["branch"]) .|> id -> (id[1],id[2])
    dm = DiscModel(
        vec(data["gen_inertia"]),
        vec(data["gen_prim_ctrl"]),
        Int64.(vec(data["gen"][:, 1])),
        findall(vec(data["bus"][:, 2]) .== 3)[1],
        coord,
        vec(data["load_freq_coef"]),
        line_list,
        -1.0 ./ data["branch"][:, 4],
        vec(data["bus"][:, 3]) / 100.0,
        vec(data["bus"][:, 9]) / 180.0 * pi,
        vec(data["gen"][:, 2]) / 100.0,
        vec(data["gen"][:, 9]) / 100.0,
        size(data["bus"], 1),
        size(data["gen"], 1),
        size(data["branch"], 1),
        )
    return dm
end
