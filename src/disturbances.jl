export add_local_disturbance!, add_local_disturbance_with_gps!


function add_local_disturbance!(model::ContModel, coord::Vector{<:Real}, dP::Real, σ::Real)::Nothing
    p = zeros(ndofs(model.dh₁))
    ch = ConstraintHandler(model.dh₁)
    db = Dirichlet(:u, Set(1:getnnodes(model.grid)), (x, _) -> dP / (σ^2 * 2 * π) * exp(-0.5 * ((x .- coord)' * (x .- coord)) / σ^2))
    add!(ch, db)
    close!(ch)
    update!(ch)
    apply!(p, ch)
    normalize_values!(p, dP, model.area, model.grid, model.dh₁, model.cellvalues)
    model.fault_nodal = p
    model.fault = (x; extrapolate=true, warn=:semi) -> interpolate(x, model.grid, model.dh₁, p, :u, extrapolate=extrapolate, warn=warn)
    return nothing
end


function add_local_disturbance_with_gps!(
    contmod::ContModel,
    gps_coord::Vector{Float64},
    dP::Float64,
    sigma::Float64,
)
    grid_coord = contmod.mesh.coord
    coord = albers_projection(reshape(gps_coord, 1, 2) ./ (180 / pi))
    coord ./= contmod.mesh.scale_factor
    add_local_disturbance!(contmod, vec(coord), dP, sigma)
    nothing
end
