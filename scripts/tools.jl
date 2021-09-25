#include("plotting.jl")
#include("ps_analysis.jl")
using JSON



function inPolygon_new(
    point::Array{Tuple{Float64, Float64}, 1},
    poly::Array{Tuple{Float64, Float64}, 1},
)
    N = length(poly)
    Np = length(point)
    b = falses(Np)
    Threads.@threads for k in 1:Np
        j = N
        for i = 1:N
            if (
                ((poly[i][2] < point[k][2]) & (poly[j][2] >= point[k][2])) |
                ((poly[j][2] < point[k][2]) & (poly[i][2] >= point[k][2]))
            )
                if (
                    poly[i][1] +
                    (point[k][2] - poly[i][2]) / (poly[j][2] - poly[i][2]) *
                    (poly[j][1] - poly[i][1]) < point[k][1]
                )
                    b[k] = !b[k]
                end
            end
            j = i
        end
    end
    return b
end



function inPolygon(p::Array{Float64,2}, poly::Array{Float64,2})
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


function alberts_projection(
    coord::Array{Float64,2};
    lon0::Float64 = 13.37616 / 180 * pi,
    lat0::Float64 = 46.94653 / 180 * pi,
    lat1::Float64 = 10 / 180 * pi,
    lat2::Float64 = 50 / 180 * pi,
    R::Float64 = 6371.0
)
    # see https://en.wikipedia.org/wiki/Albers_projection
    # R = 6.371 #Earth radius in 1000km
    n = 1 / 2 * (sin(lat1) + sin(lat2))
    theta = n * (coord[:, 1] .- lon0)
    c = cos(lat1)^2 + 2 * n * sin(lat1)
    rho = R / n * sqrt.(c .- 2 * n * sin.(coord[:, 2]))
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

    b = alberts_projection(b)
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
        println(k)
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


function back_to_2d(isgrid::BitVector, Ny::Int64, Nx::Int64, valueflat::Array{Float64,2})
    value = zeros(Ny, Nx, size(valueflat, 2))
    for t in 1:size(valueflat, 2)
    	temp = zeros(size(isgrid))
    	temp[isgrid] .= valueflat[:, t]
        value[:,:, t] = reshape(temp, Ny, Nx)
    end
    return value
end
