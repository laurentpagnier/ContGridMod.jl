export perform_dyn_sim, save_simulation

"""
    save_simulation(model::ContModel, sol::DESolution, fname::String; dt::Real=1e-2)::Nothing

Save a dynamical simulation to a .pvd file for use with Paraview.
"""
function save_simulation(model::ContModel, sol::DESolution, fname::String; dt::Real=1e-2)::Nothing
    mkpath(fname)
    pvd = paraview_collection(fname * ".pvd")
    for (i, t) in enumerate(sol.t[1]:dt:sol.t[end])
        vtk_grid(fname * "/" * fname * "-$i", model.dh₂) do vtk
            vtk_point_data(vtk, model.dh₂, sol(t))
            vtk_save(vtk)
            pvd[t] = vtk
        end
    end
    vtk_save(pvd)
    return nothing
end

"""
    perform_dyn_sim(model::ContModel, tstart::Real, tfinal::Real; solver::OrdinaryDiffEqAlgorithm=Trapezoid(), saveat::Real=0.01, abstol::Real=1e-7, reltol::Real=1e-5)::DESolution

Perform a simulation of the continuos swing equations on the model.
"""
function perform_dyn_sim(
    model::ContModel,
    tstart::Real,
    tfinal::Real;
    solver::OrdinaryDiffEqAlgorithm=Trapezoid(),
    saveat::Real=0.01,
    abstol::Real=1e-7,
    reltol::Real=1e-5
)::DESolution
    # Assemble mass matrix, stiffnes matrix, and force vector
    K = create_sparsity_pattern(model.dh₂)
    M = create_sparsity_pattern(model.dh₂)
    f = zeros(ndofs(model.dh₂))
    K, f = assemble_K!(K, f, model)
    M = assemble_M!(M, model)

    # Create intiial conditions
    ch = ConstraintHandler(model.dh₂)
    db = Dirichlet(:θ, Set(1:getnnodes(model.grid)), (x, t) -> model.θ₀(x))
    add!(ch, db)
    close!(ch)
    update!(ch, 0)
    u₀ = zeros(ndofs(model.dh₂))
    apply!(u₀, ch)

    # Create ODE problem
    function dif!(du, u, _, _)
        mul!(du, K, u)
        du .+= f
    end
    function jac!(J, _, _, _)
        J[:, :] = K
    end
    jac_sparsity = sparse(K)
    rhs = ODEFunction(dif!, mass_matrix=M, jac=jac!, jac_prototype=jac_sparsity)
    problem = ODEProblem(rhs, u₀, (tstart, tfinal))
    sol = solve(problem, solver, saveat=saveat, tstops=saveat, reltol=reltol, abstol=abstol)

    return sol
end

"""
    assemble_M!(M::SparseMatrixCSC, model::ContModel)::SparseMatrixCSC

Use the ferrite.jl assembler to create the mass matrix.
"""
function assemble_M!(M::SparseMatrixCSC, model::ContModel)::SparseMatrixCSC
    n_basefuncs_θ = getnbasefunctions(model.cellvalues)
    n_basefuncs_ω = getnbasefunctions(model.cellvalues)
    n_basefuncs = n_basefuncs_θ + n_basefuncs_ω
    θ▄, ω▄ = 1, 2

    Mₑ = PseudoBlockArray(zeros(n_basefuncs, n_basefuncs), [n_basefuncs_θ, n_basefuncs_ω], [n_basefuncs_θ, n_basefuncs_ω])

    assembler = start_assemble(M)

    for cell in CellIterator(model.dh₂)

        fill!(Mₑ, 0)


        Ferrite.reinit!(model.cellvalues, cell)

        for q_point in 1:getnquadpoints(model.cellvalues)
            x = spatial_coordinate(model.cellvalues, q_point, getcoordinates(cell))
            dΩ = getdetJdV(model.cellvalues, q_point)
            for i in 1:n_basefuncs_θ
                φᵢ = shape_value(model.cellvalues, q_point, i)
                for j in 1:n_basefuncs_θ
                    φⱼ = shape_value(model.cellvalues, q_point, j)
                    Mₑ[BlockIndex((θ▄, θ▄), (i, j))] += φᵢ ⋅ φⱼ * dΩ
                    Mₑ[BlockIndex((ω▄, ω▄), (i, j))] += model.m(x) * φᵢ ⋅ φⱼ * dΩ
                end
            end
        end

        assemble!(assembler, celldofs(cell), Mₑ)
    end
    return M
end

"""
    assemble_K!(K::SparseMatrixCSC, f::Vector{<:Real}, model::ContModel)::Tuple{SparseMatrixCSC, Vector{<:Real}}

Use the ferrite.jl assembler to create the stiffness matrix and force vector.
"""
function assemble_K!(K::SparseMatrixCSC, f::Vector{<:Real}, model::ContModel)::Tuple{SparseMatrixCSC, Vector{<:Real}}
    n_basefuncs_θ = getnbasefunctions(model.cellvalues)
    n_basefuncs_ω = getnbasefunctions(model.cellvalues)
    n_basefuncs = n_basefuncs_θ + n_basefuncs_ω
    θ▄, ω▄ = 1, 2

    Kₑ = PseudoBlockArray(zeros(n_basefuncs, n_basefuncs), [n_basefuncs_θ, n_basefuncs_ω], [n_basefuncs_θ, n_basefuncs_ω])
    fₑ = zeros(n_basefuncs)

    assembler = start_assemble(K, f)

    for cell in CellIterator(model.dh₂)
        fill!(Kₑ, 0)
        fₑ .= 0

        Ferrite.reinit!(model.cellvalues, cell)

        for q_point in 1:getnquadpoints(model.cellvalues)
            x = spatial_coordinate(model.cellvalues, q_point, getcoordinates(cell))
            b = SparseMatrixCSC(diagm([model.bx(x), model.by(x)]))
            dΩ = getdetJdV(model.cellvalues, q_point)
            for i in 1:n_basefuncs_θ
                φᵢ = shape_value(model.cellvalues, q_point, i)
                ∇φᵢ = shape_gradient(model.cellvalues, q_point, i)
                for j in 1:n_basefuncs_ω
                    φⱼ = shape_value(model.cellvalues, q_point, j)
                    ∇φⱼ = shape_gradient(model.cellvalues, q_point, j)
                    Kₑ[BlockIndex((θ▄, ω▄), (i, j))] += φᵢ ⋅ φⱼ * dΩ
                    Kₑ[BlockIndex((ω▄, θ▄), (i, j))] -= ∇φᵢ ⋅ (b * ∇φⱼ) * dΩ
                    Kₑ[BlockIndex((ω▄, ω▄), (i, j))] -= model.d(x) * φᵢ ⋅ φⱼ * dΩ
                end
                fₑ[n_basefuncs_θ+i] += (model.p(x) + model.fault(x)) * φᵢ * dΩ
            end
        end

        assemble!(assembler, celldofs(cell), Kₑ, fₑ)
    end
    return K, f
end

