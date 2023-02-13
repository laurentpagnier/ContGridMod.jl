export compute_stable_sol!, compute_stable_sol, stable_sol, stable_sol!


function stable_sol(model::ContModelFer)::Vector{Float64}
    # We need to fix the value at one grid point
    K = create_sparsity_pattern(model.dh₁)
    f = zeros(model.dh₁.ndofs.x)
    K, f = assemble_K₀!(model.K₀, model.f₀, model.cellvalues, model.dh₁, model)
    apply!(K, f, model.ch)
    model.K₀ = K
    model.f₀ = f
    θ₀ = K \ f
    return θ₀
end

function stable_sol!(model::ContModelFer)::Nothing
    θ₀ = stable_sol(model)
    model.θ₀_nodal = θ₀
    model.θ₀ = (x; extrapolate=true, warn=:semi) -> interpolate(x, model.grid, model.dh₁, θ₀, :u, extrapolate=extrapolate, warn=warn)
    return nothing
end

function assemble_K₀!(K::SparseMatrixCSC, f::Vector, cellvalues::CellScalarValues{dim}, dh::DofHandler, model::ContModelFer) where {dim}
    n_basefuncs = getnbasefunctions(cellvalues)
    Kₑ = zeros(n_basefuncs, n_basefuncs)
    fₑ = zeros(n_basefuncs)
    assembler = start_assemble(K, f)

    for cell in CellIterator(dh)
        fill!(Kₑ, 0)
        fₑ .= 0

        Ferrite.reinit!(cellvalues, cell)

        for q_point in 1:getnquadpoints(cellvalues)
            dΩ = getdetJdV(cellvalues, q_point)
            x = spatial_coordinate(cellvalues, q_point, getcoordinates(cell))
            b = SparseMatrixCSC(diagm([model.bx(x), model.by(x)]))
            for i in 1:n_basefuncs
                φᵢ = shape_value(cellvalues, q_point, i)
                ∇φᵢ = shape_gradient(cellvalues, q_point, i)
                for j in 1:n_basefuncs
                    ∇φⱼ = shape_gradient(cellvalues, q_point, j)
                    Kₑ[i, j] += ∇φᵢ ⋅ (b * ∇φⱼ) * dΩ
                end
                fₑ[i] += model.p(x) * φᵢ * dΩ
                end
        end

        assemble!(assembler, celldofs(cell), Kₑ, fₑ)
    end
    return K, f
end

function compute_stable_sol!(
    cm::ContModel,
)
    cm.th = compute_stable_sol(cm)
    nothing
end


function compute_stable_sol(
    cm::ContModel,
)
    B = sparse([cm.mesh.id_edge[:,1]; cm.mesh.id_edge[:,2]],
        [1:cm.mesh.Nedge; 1:cm.mesh.Nedge], [-ones(cm.mesh.Nedge); ones(cm.mesh.Nedge)])
    id = setdiff(1:cm.mesh.Nnode, cm.id_slack)
    I = sparse(id, 1:cm.mesh.Nnode-1, ones(cm.mesh.Nnode-1),
        cm.mesh.Nnode, cm.mesh.Nnode-1)
    Bns = B[id,:]
    Bnst = SparseMatrixCSC{Float64, Int64}(Bns')
    M = -Bns * (cm.b .* Bnst)
    return I * (M \ - cm.p[id] * cm.mesh.dx^2)
end

