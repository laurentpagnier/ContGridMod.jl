export stable_sol, stable_sol!


function stable_sol(model::ContModel)::Vector{Float64}
    K = create_sparsity_pattern(model.dh₁)
    f = zeros(model.dh₁.ndofs.x)
    assemble_K₀!(K, f, model)
    apply!(K, f, model.ch)
    return K \ f
end

function stable_sol!(model::ContModel)::Nothing
    assemble_K₀!(model.K₀, model.f₀, model)
    apply!(model.K₀, model.f₀, model.ch)
    model.θ₀_nodal = model.K₀ \ model.f₀
    model.θ₀ = (x; extrapolate=true, warn=:semi) -> interpolate(x, model.grid, model.dh₁, θ₀, :u, extrapolate=extrapolate, warn=warn)
    return nothing
end

function assemble_K₀!(K::SparseMatrixCSC, f::Vector, model::ContModel)
    n_basefuncs = getnbasefunctions(model.cellvalues)
    Kₑ = zeros(n_basefuncs, n_basefuncs)
    fₑ = zeros(n_basefuncs)
    assembler = start_assemble(K, f)

    for cell in CellIterator(model.dh₁)
        fill!(Kₑ, 0)
        fₑ .= 0

        Ferrite.reinit!(model.cellvalues, cell)

        for q_point in 1:getnquadpoints(model.cellvalues)
            dΩ = getdetJdV(model.cellvalues, q_point)
            x = spatial_coordinate(model.cellvalues, q_point, getcoordinates(cell))
            b = SparseMatrixCSC(diagm([model.bx(x), model.by(x)]))
            for i in 1:n_basefuncs
                φᵢ = shape_value(model.cellvalues, q_point, i)
                ∇φᵢ = shape_gradient(model.cellvalues, q_point, i)
                for j in 1:n_basefuncs
                    ∇φⱼ = shape_gradient(model.cellvalues, q_point, j)
                    Kₑ[i, j] += ∇φᵢ ⋅ (b * ∇φⱼ) * dΩ
                end
                fₑ[i] += model.p(x) * φᵢ * dΩ
                end
        end

        assemble!(assembler, celldofs(cell), Kₑ, fₑ)
    end
end
