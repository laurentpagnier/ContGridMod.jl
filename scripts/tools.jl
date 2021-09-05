include("plotting.jl")
include("ps_analysis.jl")


function inPolygon(p::Array{Float64,2}, poly::Array{Float64,2})
    N = size(poly, 1)
    b = Bool.(zeros(size(p, 1), 1))
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
                    b[k] = !b[k]
                end
            end
            j = i
        end
    end
    return b
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
    return b
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


function get_cont_values(
    isinside::BitMatrix,
    yrange::Array{Float64, 1},
    xrange::Array{Float64, 1},
    disc_coord::Array{Float64, 2},
    disc_values::Array{Float64, 1};
    Niter = 1000)
    
    idin = findall(isinside)
    cont_coord = zeros(length(idin), 2)
    for i = 1:length(idin)
        cont_coord[i,:] = [xrange[idin[i][2]], yrange[idin[i][1]]]
    end
    
    cont_values_vec = zeros(length(idin))
    for i = 1:size(idin, 1)
        dx = disc_coord[:,1] .- cont_coord[i,1] 
        dy = disc_coord[:,2] .- cont_coord[i,2] 
        id = argmin(dx.^2 + dy.^2)
        cont_values_vec[i] = disc_values[id]
    end
    
    cont_values = zeros(Ny, Nx)
    cont_values[idin] = cont_values_vec
    
    new_cont_values = zeros(Ny, Nx)
    tau = 0.001
    interval = 100

    @time begin
        for k in 1:Niter
            Threads.@threads for i in 2:Ny-1
                Threads.@threads for j in 2:Nx-1
                    if(isinside[i,j])
                        new_cont_values[i,j] = (1.0 - 4.0 * tau) * cont_values[i,j] + tau * (cont_values[i+1,j] +
                            cont_values[i-1,j] + cont_values[i,j+1] + cont_values[i,j-1])
                    end
                end
            end
            
            Threads.@threads for k in 1:size(n,1)
                i = Int64(n[k,1])
                j = Int64(n[k,2])
                nx = n[k,4]
                ny = n[k,3]
                if(nx == 1)
                    new_cont_values[i,j] = new_cont_values[i,j-2]
                elseif(nx == -1)
                    new_cont_values[i,j] = new_cont_values[i,j+2]
                end
                if(ny == 1)
                    new_cont_values[i,j] = new_cont_values[i-2,j]
                elseif(ny == -1)
                    new_cont_values[i,j] = new_cont_values[i+2,j]
                end
            end
            
            cont_values = copy(new_cont_values)
            
            if(mod(k,interval) == 0)
                println(k)
            end
        end
    end
    
    
    return cont_values
end


function find_time_step(isin::BitMatrix, m::Array{Float64,2}, d::Array{Float64,2}, p::Array{Float64,2},
    bx::Array{Float64,2}, by::Array{Float64,2}, dx::Float64, alpha=0.1)
    # !!!!!!! NOT FUNCTIONAL STILL SOME WORK TO DO
    Nx = size(m,2)
    Ny = size(m,1)
    bij = zeros(size(m))
    for i = 2:Ny-2
        for j = 2:Nx-2
            bij[i,j] = bx[i,j] + bx[i,j+1] + by[i-1,j] + by[i,j]
        end
    end
    gamma = d ./ m
    println(alpha*dx^2*minimum(m[isin] ./ bij[isin]))
    println(alpha*minimum(d[isin] ./ abs.(p[isin])))
end
