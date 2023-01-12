export get_params, load_discrete_model, update_model!, update_model_dm!

using HDF5
using Plots
using SparseArrays
using DifferentialEquations

function bb(dm::DiscModel, x::Union{Vector{T},Tensor{1,dim,T,dim}}, σ::Float64, bfactor::Float64)::Vector{Float64} where {T,dim}
    ds = dm.coord[dm.id_line[:, 2], :] .- dm.coord[dm.id_line[:, 1], :]
    re = zeros(2)
    for (i, row) in enumerate(eachrow(ds))
        if norm(row) == 0
            continue
        end
        b = dm.b[i] * row / norm(row) / bfactor
        t = (row[1] * (x[1] - dm.coord[dm.id_line[i, 1], 1]) + row[2] * (x[2] - dm.coord[dm.id_line[i, 1], 2])) / norm(row)^2
        dist = 0
        if 0 ≤ t ≤ 1
            dist = norm(x .- (dm.coord[dm.id_line[i, 1], :] .+ t * row))
        elseif t < 0
            dist = norm(x .- dm.coord[dm.id_line[i, 1], :])
        else
            dist = norm(x .- dm.coord[dm.id_line[i, 2], :])
        end
        re .+= abs.(b / (σ^2 * 2 * π) * exp(-0.5 * dist^2 / σ^2))
    end
    return re
end

function get_params(
    grid::Grid,
    tf::Float64,
    dm::DiscModel;
    κ::Float64=1.0,
    u_min::Float64=0.1,
    σ::Float64=0.01,
    bfactor::Float64=1.0
)::ContModelFer

    # Create the dof handler and interpolation functions
    dh₁ = DofHandler(grid)
    push!(dh₁, :u, 1)
    close!(dh₁)
    dh₂ = DofHandler(grid)
    push!(dh₂, :θ, 1)
    push!(dh₂, :ω, 1)
    close!(dh₂)
    ip = Lagrange{2,RefTetrahedron,1}()
    qr = QuadratureRule{2,RefTetrahedron}(2)
    cellvalues = CellScalarValues(qr, ip)

    # Get the area of the grid
    area = integrate(dh₁, cellvalues, (x) -> 1)

    # Create initial condition functions for diffusion process
    function d₀(x, t)
        re = 0
        for i in 1:dm.Ngen
            dif = x .- dm.coord[dm.id_gen[i], :]
            re += dm.d_gen[i] / (σ^2 * 2 * π) * exp(-0.5 * (dif' * dif) / σ^2)
        end
        for i in 1:dm.Nbus
            dif = x .- dm.coord[i, :]
            re += dm.d_load[i] / (σ^2 * 2 * π) * exp(-0.5 * (dif' * dif) / σ^2)
        end
        return max(re, u_min)
    end
    function m₀(x, t)
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
    function p₀(x, t)
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
    function bx₀(x, t)
        return bb(dm, x, σ, bfactor)[1]
    end
    function by₀(x, t)
        return bb(dm, x, σ, bfactor)[2]
    end
    d = diffusion(dh₁, cellvalues, grid, d₀, tf, κ)
    d = normalize_values!(d, sum(dm.d_load) + sum(dm.d_gen), area, grid, dh₁, cellvalues)
    m = diffusion(dh₁, cellvalues, grid, m₀, tf, κ)
    m = normalize_values!(m, sum(dm.m_gen[dm.p_gen.>0]), area, grid, dh₁, cellvalues)
    p = diffusion(dh₁, cellvalues, grid, p₀, tf, κ)
    p = normalize_values!(p, 0.0, area, grid, dh₁, cellvalues)
    bx = diffusion(dh₁, cellvalues, grid, bx₀, tf, κ)
    by = diffusion(dh₁, cellvalues, grid, by₀, tf, κ)

    return ContModelFer(
        grid,
        dh₁,
        dh₂,
        cellvalues,
        area,
        m,
        d,
        p,
        bx,
        by,
        zeros(ndofs(dh₁)), # θ₀_nodal
        zeros(ndofs(dh₁)), # fault_nodal
        (x) -> interpolate(x, grid, dh₁, m, :u),
        (x) -> interpolate(x, grid, dh₁, d, :u),
        (x) -> interpolate(x, grid, dh₁, p, :u),
        (x) -> interpolate(x, grid, dh₁, bx, :u),
        (x) -> interpolate(x, grid, dh₁, by, :u),
        (x) -> 0, # θ₀(x)
        (x) -> 0, # fault(x)
    )
end

function diffusion(dh::DofHandler, cellvalues::CellScalarValues, grid::Grid, f::Function, tf::Float64, κ::Float64)::Vector{Float64}
    M = create_sparsity_pattern(dh)
    M = assembleM!(M, dh, cellvalues)
    K = create_sparsity_pattern(dh)
    K = assembleK!(K, dh, cellvalues)

    # Create intital conditions
    ch = ConstraintHandler(dh)
    d = Dirichlet(:u, Set(1:getnnodes(grid)), f)
    add!(ch, d)
    close!(ch)
    update!(ch, 0)
    u₀ = zeros(ndofs(dh))
    apply!(u₀, ch)
    # In a future release of Ferrite this can be replaced by
    # apply_analytical!(u₀, dh, :u, f)

    function dif!(du, u, p, t)
        mul!(du, κ * K, u)
    end

    odefun = ODEFunction(dif!, mass_matrix=M, jac_prototype=K)
    problem = ODEProblem(odefun, u₀, (0.0, tf))
    sol = solve(problem, Trapezoid(), reltol=1e-5, abstol=1e-7)
    return sol.u[end]
end

function assembleK!(K::SparseMatrixCSC, dh::DofHandler, cellvalues::CellScalarValues)::SparseMatrixCSC
    n_basefuncs = getnbasefunctions(cellvalues)
    Kₑ = zeros(n_basefuncs, n_basefuncs)

    assembler = start_assemble(K)

    for cell in CellIterator(dh)
        fill!(Kₑ, 0)
        Ferrite.reinit!(cellvalues, cell)

        for q_point in 1:getnquadpoints(cellvalues)
            dΩ = getdetJdV(cellvalues, q_point)
            for i in 1:n_basefuncs
                ∇φᵢ = shape_gradient(cellvalues, q_point, i)
                for j in 1:n_basefuncs
                    ∇φⱼ = shape_gradient(cellvalues, q_point, j)
                    Kₑ[i, j] -= ∇φᵢ ⋅ ∇φⱼ * dΩ
                end
            end
        end
        assemble!(assembler, celldofs(cell), Kₑ)
    end
    return K
end

function assembleM!(M::SparseMatrixCSC, dh::DofHandler, cellvalues::CellScalarValues)::SparseMatrixCSC
    n_basefuncs = getnbasefunctions(cellvalues)
    Mₑ = zeros(n_basefuncs, n_basefuncs)

    assembler = start_assemble(M)

    for cell in CellIterator(dh)
        fill!(Mₑ, 0)
        Ferrite.reinit!(cellvalues, cell)

        for q_point in 1:getnquadpoints(cellvalues)
            dΩ = getdetJdV(cellvalues, q_point)
            for i in 1:n_basefuncs
                φᵢ = shape_value(cellvalues, q_point, i)
                for j in 1:n_basefuncs
                    φⱼ = shape_value(cellvalues, q_point, j)
                    Mₑ[i, j] += φᵢ ⋅ φⱼ * dΩ
                end
            end
        end
        assemble!(assembler, celldofs(cell), Mₑ)
    end
    return M
end

function integrate(dh::DofHandler, cellvalues::CellScalarValues, f::Function)::Float64
    n_basefuncs = getnbasefunctions(cellvalues)
    int = 0.0
    for cell in CellIterator(dh)
        Ferrite.reinit!(cellvalues, cell)
        for q_point in 1:getnquadpoints(cellvalues)
            x = spatial_coordinate(cellvalues, q_point, getcoordinates(cell))
            dΩ = getdetJdV(cellvalues, q_point)
            for i in 1:n_basefuncs
                φᵢ = shape_value(cellvalues, q_point, i)
                int += f(x) * φᵢ * dΩ
            end
        end
    end
    return int
end

function interpolate(x::Tensor{1,dim,T,dim}, grid::Grid, dh::DofHandler, u::Vector{Float64}, fname::Symbol; off::Float64=0.0)::Float64 where {dim,T}
    ph = PointEvalHandler(grid, [x])
    return get_point_values(ph, dh, u, fname)[1] + off
end

function normalize_values!(u::Vector{Float64}, value::Float64, area::Float64, grid::Grid, dh::DofHandler, cellvalues::CellScalarValues)::Vector{Float64}
    utot = integrate(dh, cellvalues, x -> interpolate(x, grid, dh, u, :u))
    ch = ConstraintHandler(dh)
    d = Dirichlet(:u, Set(1:getnnodes(grid)), (x, t) -> interpolate(x, grid, dh, u, :u, off=(value - utot) / area))
    add!(ch, d)
    close!(ch)
    update!(ch, 0)
    apply!(u, ch)
    return u
end

function uniform_param(
    mesh::Mesh,
    scale_factor::Float64,
    dataname::String,
)

    dm = load_discrete_model(dataname, scale_factor)

    N = size(mesh.coord, 1)

    inertia = sum(dm.m_gen[dm.p_gen.>0])
    damping = sum(dm.d_load) + sum(dm.d_gen[dm.p_gen.>0])

    m = inertia / (N * mesh.dx^2)
    d = damping / (N * mesh.dx^2)
    b = 50
    pl = (sum(dm.load) / sum(pl) / mesh.dx^2) .* pl
    pg = (sum(dm.gen) / sum(pg) / mesh.dx^2) .* pg

    # ensure that the integral of the power injection is 0, i.e
    # that generation matches the demand.

    temp = sparse([mesh.inc_mat[:, 1]; mesh.inc_mat[:, 2]], [1:n; 1:n], 0.5 * [ones(n); ones(n)])
    b = temp' * b
    temp = sparse([mesh.inc_mat[:, 1]; mesh.inc_mat[:, 2]], [1:n; 1:n], [-ones(n); ones(n)])
    xi = -temp * (b .* temp')
    minv = m .^ (-1)
    gamma = d .* minv

    id_slack = find_node(mesh, reshape(dm.coord[dm.id_slack, :], 1, 2))[1]

    p = pg - pl
    p = p .- sum(p) / length(p)
    dp = zeros(size(p))
    return cm = ContModel(
        mesh,
        id_slack,
        minv,
        gamma,
        p,
        dp,
        xi,
        b,
        m,
        d,
        zeros(size(d)),
        scale_factor,
        dmax,
        Niter,
        tau,
        patch,
        min_factor
    )
end


function get_params(
    mesh::Mesh,
    dataname::String;
    dmax::Float64=100.0,
    Niter::Int64=5000,
    tau::Float64=0.001,
    patch::Float64=0.001,
    min_factor::Float64=0.1
)
    # this function uses the heat equation to "diffuse" the discrete
    # parameters

    # load the discrete model
    dm = load_discrete_model(dataname, mesh.scale_factor)
    #idin = findall(mesh.isinside)
    N = size(mesh.coord, 1)
    m = zeros(N)
    d = zeros(N)
    pl = zeros(N)
    pg = zeros(N)
    b = zeros(N)

    Threads.@threads for g in 1:dm.Ngen
        if (dm.p_gen[g] > 0.0)
            k = argmin((mesh.coord[:, 1] .- dm.coord[dm.id_gen[g], 1]) .^ 2 +
                       (mesh.coord[:, 2] .- dm.coord[dm.id_gen[g], 2]) .^ 2)
            m[k] += dm.m_gen[g]
            d[k] += dm.d_gen[g]
            pg[k] += dm.p_gen[g]
        end
    end

    Threads.@threads for l in 1:dm.Nbus
        k = argmin((mesh.coord[:, 1] .- dm.coord[l, 1]) .^ 2 +
                   (mesh.coord[:, 2] .- dm.coord[l, 2]) .^ 2)
        d[k] += dm.d_load[l]
        pl[k] += dm.p_load[l]
    end

    bm = sum(dm.b) / length(dm.b)
    Threads.@threads for l in 1:dm.Nline
        x2 = dm.coord[dm.id_line[l, 2], 2]
        x1 = dm.coord[dm.id_line[l, 1], 2]
        y2 = dm.coord[dm.id_line[l, 2], 1]
        y1 = dm.coord[dm.id_line[l, 1], 1]
        dx_l = x2 - x1
        dy_l = y2 - y1
        ds2 = (dy_l^2 + dx_l^2)
        if (dx_l != 0 && dy_l != 0) # if not a transformer
            y = mesh.coord[:, 1] # it's "less precise" but way simpler to implement 
            x = mesh.coord[:, 2] # using node locations instead of centers of lines, 
            beta = (dx_l .* (y1 .- y) + dy_l .* (x .- x1)) ./ ds2
            alpha = (dx_l .* (x .- x1) + dy_l .* (y .- y1)) ./ ds2
            in_seg = (0 .< alpha) .& (alpha .< 1)
            dist = abs.(beta) .* in_seg .* sqrt(ds2) + # if close to the segment
                   .!in_seg .* min.(sqrt.((x .- x1) .^ 2 + (y .- y1) .^ 2), # if close to the ends
                sqrt.((x .- x2) .^ 2 + (y .- y2) .^ 2))
            if (dm.b[l] < 2 * bm) # if b is too large discard it
                b[dist.<dmax] .+= dm.b[l] * mesh.dx^2 * patch
            end
        end
    end

    #n = size(mesh.id_line, 1)
    N = size(mesh.coord, 1)

    m = heat_diff(mesh, m, Niter=Niter, tau=tau)
    d = heat_diff(mesh, d, Niter=Niter, tau=tau)
    pl = heat_diff(mesh, pl, Niter=Niter, tau=tau)
    pg = heat_diff(mesh, pg, Niter=Niter, tau=tau)
    b = heat_diff(mesh, b, Niter=Niter, tau=tau)

    # asign minimal values to the quantities
    m .= max.(m, min_factor * sum(m) / length(m))
    d .= max.(d, min_factor * sum(d) / length(m))

    # ensure that the integrals of the different quantities over
    # the medium is equivalent to their sum over the discrete model
    m = (sum(dm.m_gen[dm.p_gen.>0]) / sum(m) / mesh.dx^2) .* m
    d = ((sum(dm.d_gen[dm.p_gen.>0]) + sum(dm.d_load)) / sum(d) / mesh.dx^2) .* d
    pl = (sum(dm.p_load) / sum(pl) / mesh.dx^2) .* pl
    pg = (sum(dm.p_gen) / sum(pg) / mesh.dx^2) .* pg

    # ensure that the integral of the power injection is 0, i.e
    # that generation matches the demand.

    temp = sparse([mesh.id_edge[:, 1]; mesh.id_edge[:, 2]],
        [1:mesh.Nedge; 1:mesh.Nedge], 0.5 * [ones(mesh.Nedge); ones(mesh.Nedge)])
    b = temp' * b

    id_slack = find_node(mesh, reshape(dm.coord[dm.id_slack, :], 1, 2))[1]

    p = pg - pl
    p = p .- sum(p) / length(p)
    dp = zeros(size(p))
    return cm = ContModel(
        mesh,
        id_slack,
        p,
        dp,
        b,
        m,
        d,
        zeros(size(d)), # th
    )
end


function heat_diff(
    mesh::Mesh,
    q::Array{Float64,1};
    Niter::Int64=5000,
    tau::Float64=0.001
)
    # parameters over the lattice with tau=kappa*dt
    #n = size(mesh.inc_mat, 1)
    #N = size(mesh.coord, 1)
    # incidence matrix
    temp = sparse([mesh.id_edge[:, 1]; mesh.id_edge[:, 2]],
        [1:mesh.Nedge; 1:mesh.Nedge], [-ones(mesh.Nedge); ones(mesh.Nedge)])
    # - Laplacian
    xi = -temp * temp'

    # Crank-Nicolson
    I = sparse(1:mesh.Nnode, 1:mesh.Nnode, ones(mesh.Nnode))
    A = 2.0 .* I - tau .* xi / mesh.dx^2
    B = 2.0 .* I + tau .* xi / mesh.dx^2
    for _ in 1:Niter
        q = A \ (B * q) #way slower when dx -> 0
        #x = A \ x
        #gmres!(x, A , B * x)
    end

    return q
end


function load_discrete_model(
    dataname::String,
    scaling_factor::Float64
)
    data = h5read(dataname, "/")
    coord = albers_projection(data["bus_coord"] ./ (180 / pi))
    coord ./= scaling_factor
    dm = DiscModel(
        vec(data["gen_inertia"]),
        vec(data["gen_prim_ctrl"]),
        Int64.(vec(data["gen"][:, 1])),
        findall(vec(data["bus"][:, 2]) .== 3)[1],
        coord,
        vec(data["load_freq_coef"]),
        Int64.(data["branch"][:, 1:2]),
        1.0 ./ data["branch"][:, 4],
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


function update_model!(
    contmod::ContModel,
    filename::String,
)
    # Not foolproof, expect the data to have the right size
    data = h5read(filename, "/")
    data_keys = collect(keys(data))
    model_keys = String.(collect(fieldnames(ContGridMod.ContModel)))
    for k in data_keys
        if k in model_keys
            setfield!(contmod, Symbol(k), data[k])
        end
    end
end

function update_model_dm!(
    contmod::ContModel,
    dm::DiscModel;
    Niter::Int64=5000,
    tau::Float64=0.001,
    min_factor::Float64=0.1
)
    # Update the dynamical parameters of the model from a discrete model
    N = size(contmod.mesh.coord, 1)

    m = zeros(N)
    d = zeros(N)
    pl = zeros(N)
    pg = zeros(N)

    Threads.@threads for g in 1:dm.Ngen
        if (dm.p_gen[g] > 0.0)
            k = argmin((contmod.mesh.coord[:, 1] .- dm.coord[dm.id_gen[g], 1]) .^ 2 +
                       (contmod.mesh.coord[:, 2] .- dm.coord[dm.id_gen[g], 2]) .^ 2)
            m[k] += dm.m_gen[g]
            d[k] += dm.d_gen[g]
            pg[k] += dm.p_gen[g]
        end
    end

    Threads.@threads for l in 1:dm.Nbus
        k = argmin((contmod.mesh.coord[:, 1] .- dm.coord[l, 1]) .^ 2 +
                   (contmod.mesh.coord[:, 2] .- dm.coord[l, 2]) .^ 2)
        d[k] += dm.d_load[l]
        pl[k] += dm.p_load[l]
    end
    @time begin
        m = heat_diff(contmod.mesh, m, Niter=Niter, tau=tau)
        d = heat_diff(contmod.mesh, d, Niter=Niter, tau=tau)
        pl = heat_diff(contmod.mesh, pl, Niter=Niter, tau=tau)
        pg = heat_diff(contmod.mesh, pg, Niter=Niter, tau=tau)
    end

    # asign minimal values to the quantities
    m .= max.(m, min_factor * sum(m) / length(m))
    d .= max.(d, min_factor * sum(d) / length(m))

    # ensure that the integrals of the different quantities over
    # the medium is equivalent to their sum over the discrete model
    m = (sum(dm.m_gen[dm.p_gen.>0]) / sum(m) / contmod.mesh.dx^2) .* m
    d = ((sum(dm.d_gen[dm.p_gen.>0]) + sum(dm.d_load)) / sum(d) / contmod.mesh.dx^2) .* d
    pl = (sum(dm.p_load) / sum(pl) / contmod.mesh.dx^2) .* pl
    pg = (sum(dm.p_gen) / sum(pg) / contmod.mesh.dx^2) .* pg

    # ensure that the integral of the power injection is 0, i.e
    # that generation matches the demand.
    p = pg - pl
    p = p .- sum(p) / length(p)

    # Update the paramters
    contmod.m = m
    contmod.d = d
    contmod.p = p
end
