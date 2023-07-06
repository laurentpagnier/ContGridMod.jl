export update_model!, albers_projection, import_border, to_jld2, from_jld2

"""
    update_model!(model::ContModel, u_name::Symbol, u::Vector{<:Real})::Nothing

Update the continuos model from a vector of nodal values. It is sufficient to pass the symbol name without the _nodal suffix.
"""
function update_model!(model::ContModel, u_name::Symbol, u::Vector{<:Real})::Nothing
    if size(u) != (getnnodes(model.grid),)
        throw(DimensionMismatch("The size of the vector does not match the number of nodes."))
    end
    setproperty!(model, Symbol(string(u_name) * "_nodal"), u)
    setproperty!(model, u_name, (x; extrapolate=true, warn=:semi) -> interpolate(x, model.grid, model.dh₁, u, :u, extrapolate=extrapolate, warn=warn))
    return nothing
end

"""
    update_model(model::ContModel,u_name::Symbol, dm::DiscModel, tf::Real[, κ::Real=1.0, u_min::Real=0.1, σ::Real=0.01, bfactor::Real=1.0, bmin::Real=1000.0])::Nothing

Update the continuous model from a discrete model. The paramter are distributed using a diffusion procees (similar to the get_params function)
"""
function update_model!(
    model::ContModel,
    u_name::Symbol,
    dm::DiscModel,
    tf::Real;
    κ::Real=1.0,
    u_min::Real=0.1,
    σ::Real=0.01,
    bfactor::Real=1.0,
    bmin::Real=1000.0
)::Nothing

    ## This whole function and the get_params should be refactored ##
    area = integrate(model.dh₁, model.cellvalues, (x) -> 1)
    if u_name == :d
        function d₀(x, _)
            re = 0
            for i in 1:dm.Ngen
                if dm.p_gen == 0
                    continue
                end
                dif = x .- dm.coord[dm.id_gen[i], :]
                re += dm.d_gen[i] / (σ^2 * 2 * π) * exp(-0.5 * (dif' * dif) / σ^2)
            end
            for i in 1:dm.Nbus
                dif = x .- dm.coord[i, :]
                re += dm.d_load[i] / (σ^2 * 2 * π) * exp(-0.5 * (dif' * dif) / σ^2)
            end
            return max(re, u_min)
        end
        d = diffusion(model.dh₁, model.cellvalues, model.grid, d₀, tf, κ)
        d = normalize_values!(d, sum(dm.d_load) + sum(dm.d_gen[dm.p_gen.>0]), area, model.grid, model.dh₁, model.cellvalues)
        update_model!(model, u_name, d)
    elseif u_name == :m
        function m₀(x, _)
            re = 0
            for i in 1:dm.Ngen
                if dm.p_gen[i] == 0
                    continue
                end
                dif = x .- dm.coord[dm.id_gen[i], :]
                re += dm.m_gen[i] / (σ^2 * 2 * π) * exp(-0.5 * (dif' * dif) / σ^2)
            end
            return max(re, u_min)
        end
        m = diffusion(model.dh₁, model.cellvalues, model.grid, m₀, tf, κ)
        m = normalize_values!(m, sum(dm.m_gen[dm.p_gen.>0]), area, model.grid, model.dh₁, model.cellvalues)
        update_model!(model, u_name, m)
    elseif u_name == :p
        function p₀(x, _)
            re = 0
            for i in 1:dm.Ngen
                dif = x .- dm.coord[dm.id_gen[i], :]
                re += dm.p_gen[i] / (σ^2 * 2 * π) * exp(-0.5 * (dif' * dif) / σ^2)
            end
            for i in 1:dm.Nbus
                dif = x .- dm.coord[i, :]
                re -= dm.p_load[i] / (σ^2 * 2 * π) * exp(-0.5 * (dif' * dif) / σ^2)
            end
            return re
        end
        p = diffusion(model.dh₁, model.cellvalues, model.grid, p₀, tf, κ)
        p = normalize_values!(p, 0.0, area, model.grid, model.dh₁, model.cellvalues, mode="off")
        update_model!(model, u_name, p)
    elseif u_name == :bx
        function bx₀(x, _)
            return bb(dm, x, σ, bfactor)[1]
        end
        bx = diffusion(model.dh₁, model.cellvalues, model.grid, bx₀, tf, κ)
        bx .= max.(bx, bmin)
        update_model!(model, u_name, bx)
    elseif u_name == :by
        function by₀(x, _)
            return bb(dm, x, σ, bfactor)[2]
        end
        by = diffusion(model.dh₁, model.cellvalues, model.grid, by₀, tf, κ)
        by .= max.(by, bmin)
        update_model!(model, u_name, by)
    end

    return nothing
end

"""
    albers_projection(coord::Array{<:Real,2}; # as latitude, longitude  lon0::Real=13.37616 / 180 * pi, lat0::Real=46.94653 / 180 * pi, lat1::Real=10 / 180 * pi, lat2::Real=50 / 180 * pi, R::Real=6371.0)::Array{<:Real,2}

Apply the Albers projection to a vector of coordinates. The coordinates need to be given as latitude, longitude.
See https://en.wikipedia.org/wiki/Albers_projection
"""
function albers_projection(
    coord::Array{<:Real,2};
    lon0::Real=13.37616 / 180 * pi,
    lat0::Real=46.94653 / 180 * pi,
    lat1::Real=10 / 180 * pi,
    lat2::Real=50 / 180 * pi,
    R::Real=6371.0
)::Array{<:Real,2}
    n = 1 / 2 * (sin(lat1) + sin(lat2))
    theta = n * (coord[:, 2] .- lon0)
    c = cos(lat1)^2 + 2 * n * sin(lat1)
    rho = R / n * sqrt.(c .- 2 * n * sin.(coord[:, 1]))
    rho0 = R / n * sqrt.(c - 2 * n * sin(lat0))
    x = rho .* sin.(theta)
    y = rho0 .- rho .* cos.(theta)
    return [x y]
end

"""
    import_border(filename::String)::Tuple{Array{2,<:Real},Real}

Import border from a json file, apply the Albers projection and rescale it such that the longest dimension is 1.
The coordinates need to be given as latitude, longitude.
"""
function import_border(
    filename::String
)::Tuple{Array{<:Real,2},Real}
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
        maximum(b[:, 1]) - minimum(b[:, 1]),
        maximum(b[:, 2]) - minimum(b[:, 2])
    )
    return b / scale_factor, scale_factor
end

"""
    to_jld2(fn::String, model::Union{ContModel,DiscModel})::nothing

Save a continuous or discrete model to a jld2 file.
"""
function to_jld2(
    fn::String,
    model::Union{ContModel,DiscModel}
)::nothing
    tmp = Dict(string(key) => getfield(model, key) for key ∈ fieldnames(typeof(model)))
    tmp["type"] = typeof(model)
    if (typeof(model) == ContModel)
        delete!(tmp, "m")
        delete!(tmp, "d")
        delete!(tmp, "p")
        delete!(tmp, "bx")
        delete!(tmp, "by")
        delete!(tmp, "θ₀")
        delete!(tmp, "fault")
    end
    save(fn, tmp)
end

"""
    from_jld2(fn::String)::Union{ContModel, DiscModel}

Load a continuous or discrete model from a jld2 file.
"""
function from_jld2(fn::String)::Union{ContModel,DiscModel}
    tmp = load(fn)
    if (tmp["type"] == ContModel)
        return ContModel(
            tmp["grid"],
            tmp["dh₁"],
            tmp["dh₂"],
            tmp["cellvalues"],
            tmp["area"],
            tmp["m_nodal"],
            tmp["d_nodal"],
            tmp["p_nodal"],
            tmp["bx_nodal"],
            tmp["by_nodal"],
            tmp["θ₀_nodal"],
            tmp["fault_nodal"],
            (x; extrapolate=true, warn=:semi) -> interpolate(x, tmp["grid"], tmp["dh₁"], tmp["m_nodal"], :u, extrapolate=extrapolate, warn=warn),
            (x; extrapolate=true, warn=:semi) -> interpolate(x, tmp["grid"], tmp["dh₁"], tmp["d_nodal"], :u, extrapolate=extrapolate, warn=warn),
            (x; extrapolate=true, warn=:semi) -> interpolate(x, tmp["grid"], tmp["dh₁"], tmp["p_nodal"], :u, extrapolate=extrapolate, warn=warn),
            (x; extrapolate=true, warn=:semi) -> interpolate(x, tmp["grid"], tmp["dh₁"], tmp["bx_nodal"], :u, extrapolate=extrapolate, warn=warn),
            (x; extrapolate=true, warn=:semi) -> interpolate(x, tmp["grid"], tmp["dh₁"], tmp["by_nodal"], :u, extrapolate=extrapolate, warn=warn),
            (x; extrapolate=true, warn=:semi) -> interpolate(x, tmp["grid"], tmp["dh₁"], tmp["θ₀_nodal"], :u, extrapolate=extrapolate, warn=warn),
            (x; extrapolate=true, warn=:semi) -> interpolate(x, tmp["grid"], tmp["dh₁"], tmp["fault_nodal"], :u, extrapolate=extrapolate, warn=warn),
            tmp["ch"],
            tmp["K₀"],
            tmp["f₀"]
        )
    elseif (tmp["type"] == DiscModel)
        return DiscModel(
            tmp["m_gen"],
            tmp["d_gen"],
            tmp["id_gen"],
            tmp["id_slack"],
            tmp["coord"],
            tmp["d_load"],
            tmp["id_line"],
            tmp["b"],
            tmp["p_load"],
            tmp["th"],
            tmp["p_gen"],
            tmp["max_gen"],
            tmp["Nbus"],
            tmp["Ngen"],
            tmp["Nline"],
        )
    else
        throw(ArgumentError("The provided file does not have a type entry matching ContGridMod.ContModelFer or ContGridMod.DiscModel"))
    end
end
