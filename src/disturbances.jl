export add_local_disturbance!, add_local_disturbance_with_gps!

function add_local_disturbance!(model::ContModelFer, coord::Vector{Float64}, dP::Real, σ::Real)::Nothing
    p = zeros(ndofs(model.dh₁))
    ch = ConstraintHandler(model.dh₁)
    db = Dirichlet(:u, Set(1:getnnodes(model.grid)), (x, t) -> exponential2D(x, coord, dP, σ))
    add!(ch, db)
    close!(ch)
    update!(ch)
    apply!(p, ch)
    normalize_values!(p, dP, model.area, model.grid, model.dh₁, model.cellvalues)
    model.fault_nodal = p
    model.fault = (x; extrapolate=true, warn=:semi) -> interpolate(x, model.grid, model.dh₁, p, :u, extrapolate=extrapolate, warn=warn)
    return nothing
end

function add_local_disturbance!(
    contmod::ContModel,
    coord::Vector{Float64},
    dP::Float64,
    sigma::Float64,
)
    grid_coord = contmod.mesh.coord
    dp = zeros(size(grid_coord,1))
    for k = 1:size(grid_coord,1)
        dp[k] = exp(-((grid_coord[k,1] - coord[1])^2 + (grid_coord[k,2] - coord[2])^2) / 2 / sigma^2)
    end
    contmod.dp = dP .* dp ./ sum(dp) / contmod.mesh.dx^2
    nothing
end


function add_local_disturbance_with_gps!(
    contmod::ContModel,
    gps_coord::Vector{Float64},
    dP::Float64,
    sigma::Float64,
)
    grid_coord = contmod.mesh.coord
    coord = albers_projection(reshape(gps_coord,1,2) ./ (180 / pi) )
    coord ./= contmod.mesh.scale_factor
    add_local_disturbance!(contmod, vec(coord), dP, sigma)
    nothing
end
