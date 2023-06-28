export perform_dyn_sim, perform_dyn_sim_backward_euler, perform_dyn_sim_crank_nicolson, perform_dyn_sim_forward, perform_dyn_sim_runge_kutta, save_simulation

function save_simulation(model::ContModelFer, sol::DiffEqBase.DESolution, fname::String)::Nothing
    mkpath(fname)
    pvd = paraview_collection(fname * ".pvd")
    for (i, t) in enumerate(sol.t)
        vtk_grid(fname * "/" * fname * "-$i", model.dh₂) do vtk
            vtk_point_data(vtk, model.dh₂, sol.u[i])
            vtk_save(vtk)
            pvd[t] = vtk
        end
    end
    vtk_save(pvd)
    return nothing
end

function perform_dyn_sim(
    model::ContModelFer,
    tf::Float64;
    solver::DiffEqBase.DEAlgorithm=Trapezoid(),
    saveat::Real=0.01,
    abstol::Float64=1e-7,
    reltol::Float64=1e-5
)::DiffEqBase.DESolution
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
    jac_sparsity = sparse(K)
    rhs = ODEFunction(dif!, mass_matrix=M, jac_prototype=jac_sparsity)
    problem = ODEProblem(rhs, u₀, (0, tf), (K, f))
    sol = solve(problem, solver, saveat=saveat, tstops=saveat, reltol=reltol, abstol=abstol)

    return sol
end

function dif!(du, u, p, t)
    mul!(du, p[1], u)
    du .+= p[2]
end

function assemble_M!(M::SparseMatrixCSC, model::ContModelFer)

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

function assemble_K!(K::SparseMatrixCSC, f::Vector, model::ContModelFer)
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


function perform_dyn_sim(
    contmod::ContModel;
    interval::Int64,
    Ndt::Int64,
    dt::Float64,
    method::String="crank-nicolson"
)

    if (method == "crank-nicolson" || method == "cn")
        return ts, thetas, omegas = perform_dyn_sim_crank_nicolson(
            contmod; interval, Ndt, dt)
    elseif (method == "backward-euler" || method == "be")
        return ts, thetas, omegas = perform_dyn_sim_backward_euler(
            contmod; interval, Ndt, dt)
    elseif (method == "runge-kutta" || method == "rk")
        return ts, thetas, omegas = perform_dyn_sim_runge_kutta(
            contmod; interval, Ndt, dt)
    elseif (method == "forward" || method == "cf")
        return ts, thetas, omegas = perform_dyn_sim_forward(
            contmod; interval, Ndt, dt)
    else
        println("Method not found, '", method, "' is not implemented.")
    end

end


function perform_dyn_sim_forward(
    contmod::ContModel;
    interval::Int64=100,
    Ndt::Int64=15000,
    dt::Float64=0.0001
)
    println("Total time: ", dt * Ndt)
    N = sum(contmod.isgrid)
    M = 1 + Int64(ceil(Ndt / interval))
    omegas = zeros(N, M)
    thetas = zeros(N, M)
    th_new = zeros(N)
    ts = zeros(M)

    th_old = copy(contmod.th[contmod.isgrid])
    th = copy(contmod.th[contmod.isgrid])
    omegas[:, 1] = zeros(size(th))
    thetas[:, 1] = copy(contmod.th[contmod.isgrid])

    minv = 1 ./ contmod.m
    gamma = contmod.d ./ contmod.m
    Nedge = contmod.mesh.Nedge
    temp = sparse([contmod.mesh.id_edge[:, 1]; contmod.mesh.id_edge[:, 2]],
        [1:Nedge; 1:Nedge], [-ones(Nedge); ones(Nedge)])
    L = -temp * (contmod.b .* temp')

    chi = 1 ./ (1 .+ gamma * dt / 2)

    A = (2 * chi - dt^2 / dx^2 * minv)
    B = (1 .- gamma * dt / 2) .* chi
    C = dt^2 * chi .* minv .* (contmod.p + contmod.dp)

    @time begin
        for t in 1:Ndt
            th_new = 2 * chi .* th +
                     dt^2 / dx^2 .* minv .* chi .* (L * th) -
                     B .* th_old + C
            if (mod(t, interval) == 0)
                omegas[:, Int64(t / interval)+1] = (th_new - th) / dt
                thetas[:, Int64(t / interval)+1] = th_new
                ts[Int64(t / interval)+1] = t * dt
                pprintln("NIter: $t Avg. Omega: $(sum(omegas[:, Int64(t/interval) + 1]) / N)")
            end
            th_old = copy(th)
            th = copy(th_new)
        end
    end

    return ts, thetas, omegas
end


function DP54(M::SparseMatrixCSC{Float64,Int64}, x::Vector{Float64}, b::Vector{Float64}; dt=0.001)
    # Dormand–Prince method
    a21 = 1.0 / 5.0
    a31 = 3.0 / 40.0
    a32 = 9.0 / 40.0
    a41 = 44.0 / 45.0
    a42 = -56.0 / 15.0
    a43 = 32.0 / 9.0
    a51 = 19372.0 / 6561.0
    a52 = -25360.0 / 2187.0
    a53 = 64448.0 / 6561.0
    a54 = -212 / 729.0
    a61 = 9017.0 / 3168.0
    a62 = -355.0 / 33.0
    a63 = 46732.0 / 5247.0
    a64 = 49.0 / 176.0
    a65 = -5103.0 / 18656.0
    a71 = 35.0 / 384.0
    a72 = 0.0
    a73 = 500.0 / 1113.0
    a74 = 125.0 / 192.0
    a75 = -2187.0 / 6784.0
    a76 = 11.0 / 84.0
    b11 = 35.0 / 384
    b12 = 0.0
    b13 = 500.0 / 1113.0
    b14 = 125.0 / 192.0
    b15 = -2187.0 / 6784.0
    b16 = 11.0 / 84.0
    b17 = 0.0
    b21 = 5179.0 / 57600.0
    b22 = 0.0
    b23 = 7571.0 / 16695.0
    b24 = 393.0 / 640.0
    b25 = -92097 / 339200.0
    b26 = 187.0 / 2100.0
    b27 = 1.0 / 40.0

    k1 = M * x + b
    k2 = M * (x + dt * a21 * k1) + b
    k3 = M * (x + dt * (a31 * k1 + a32 * k2)) + b
    k4 = M * (x + dt * (a41 * k1 + a42 * k2 + a43 * k3)) + b
    k5 = M * (x + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)) + b
    k6 = M * (x + dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5)) + b
    k7 = M * (x + dt * (a71 * k1 + a72 * k2 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6)) + b
    y1 = x + dt * (b11 * k1 + b12 * k2 + b13 * k3 + b14 * k4 + b15 * k5 + b16 * k6 + b17 * k7)
    y2 = x + dt * (b21 * k1 + b22 * k2 + b23 * k3 + b24 * k4 + b25 * k5 + b26 * k6 + b27 * k7)
    return y1, y1 - y2
end


function RK4(M::SparseMatrixCSC{Float64,Int64}, x::Vector{Float64}, b::Vector{Float64}; dt=0.001)
    k1 = M * x + b
    k2 = M * (x + dt * k1 / 2.0) + b
    k3 = M * (x + dt * k2 / 2.0) + b
    k4 = M * (x + dt * k3) + b
    y = x + dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
end


function perform_dyn_sim_runge_kutta(
    contmod::ContModel;
    interval::Int64=10,
    Ndt::Int64=1000,
    dt::Float64=0.05
)
    println("Total time: ", dt * Ndt)
    N = contmod.mesh.Nnode
    M = 1 + Int64(ceil(Ndt / interval))
    omegas = zeros(N, M)
    thetas = zeros(N, M)
    ts = zeros(M)

    x = [copy(contmod.th); zeros(N)]
    omegas[:, 1] = zeros(N)
    thetas[:, 1] = copy(contmod.th)

    ts = zeros(1 + Int64(ceil(Ndt / interval)))

    Nnode = contmod.mesh.Nnode
    Nedge = contmod.mesh.Nedge
    temp = sparse([contmod.mesh.id_edge[:, 1]; contmod.mesh.id_edge[:, 2]],
        [1:Nedge; 1:Nedge], [-ones(Nedge); ones(Nedge)])
    L = temp * (contmod.b .* temp') / contmod.mesh.dx^2

    minv = 1 ./ contmod.m
    G = sparse(1:Nnode, 1:Nnode, contmod.d .* minv)

    I = sparse(1:Nnode, 1:Nnode, ones(Nnode))

    A = [spzeros(Nnode, Nnode) I
        -minv.*L -G]

    rhs = [zeros(Nnode); minv .* (contmod.p + contmod.dp)]

    @time begin
        for t in 1:Ndt
            x = RK4(A, x, rhs, dt=dt)
            if mod(t, interval) == 0
                thetas[:, Int64(t / interval)+1] = x[1:N]
                omegas[:, Int64(t / interval)+1] = x[N+1:end]
                ts[Int64(t / interval)+1] = t * dt
                println("NIter: $t Avg. Omega: $(sum(omegas[:, Int64(t/interval) + 1]) / N)")
            end
        end
    end
    return ts, thetas, omegas
end


function perform_dyn_sim_crank_nicolson(
    contmod::ContModel;
    interval::Int64=10,
    Ndt::Int64=1000,
    dt::Float64=0.05
)
    println("Total time: ", dt * Ndt)
    Nnode = contmod.mesh.Nnode
    Ninterval = 1 + Int64(ceil(Ndt / interval))
    omegas = zeros(Nnode, Ninterval)
    thetas = zeros(Nnode, Ninterval)
    ts = zeros(Ninterval)

    x = [copy(contmod.th); zeros(Nnode)]
    omegas[:, 1] = zeros(Nnode)
    thetas[:, 1] = copy(contmod.th)

    ts = zeros(1 + Int64(ceil(Ndt / interval)))

    minv = 1 ./ contmod.m
    gamma = contmod.d ./ contmod.m
    Nedge = contmod.mesh.Nedge
    inc_mat = sparse([contmod.mesh.id_edge[:, 1]; contmod.mesh.id_edge[:, 2]],
        [1:Nedge; 1:Nedge], [-ones(Nedge); ones(Nedge)])
    L = inc_mat * (contmod.b .* inc_mat') / contmod.mesh.dx^2
    G = sparse(1:Nnode, 1:Nnode, gamma)

    I = sparse(1:Nnode, 1:Nnode, ones(Nnode))
    A = [I -dt/2*I
        dt/2*minv.*L (I+dt/2*G)]
    B = [I dt/2*I
        -dt/2*minv.*L (I-dt/2*G)]
    C = [zeros(Nnode); dt * minv .* (contmod.p + contmod.dp)]

    @time begin
        for t in 1:Ndt
            x = A \ (B * x + C) # way slower when dx -> 0
            #gmres!(x, A , B * x + C)
            if mod(t, interval) == 0
                thetas[:, Int64(t / interval)+1] = x[1:Nnode]
                omegas[:, Int64(t / interval)+1] = x[Nnode+1:end]
                ts[Int64(t / interval)+1] = t * dt
                println("NIter: $t Avg. Omega: $(sum(omegas[:, Int64(t/interval) + 1]) / Nnode)")
            end
        end
    end
    return ts, thetas, omegas
end


function perform_dyn_sim_backward_euler(
    contmod::ContModel;
    interval::Int64=10,
    Ndt::Int64=1000,
    dt::Float64=0.05
)
    println("Total time: ", dt * Ndt)
    N = contmod.mesh.Nnode
    M = 1 + Int64(ceil(Ndt / interval))
    omegas = zeros(N, M)
    thetas = zeros(N, M)
    ts = zeros(M)

    x = [copy(contmod.th); zeros(N)]
    omegas[:, 1] = zeros(N)
    thetas[:, 1] = copy(contmod.th)

    ts = zeros(1 + Int64(ceil(Ndt / interval)))

    I = sparse(1:N, 1:N, ones(N))
    minv = 1 ./ contmod.m
    gamma = contmod.d ./ contmod.m
    Nedge = contmod.mesh.Nedge
    inc_mat = sparse([contmod.mesh.id_edge[:, 1]; contmod.mesh.id_edge[:, 2]],
        [1:Nedge; 1:Nedge], [-ones(Nedge); ones(Nedge)])
    L = inc_mat * (contmod.b .* inc_mat') / contmod.mesh.dx^2

    A = [I -dt*I
        dt*minv.*L (I+dt*sparse(1:N, 1:N, gamma))]
    B = [zeros(N); dt * minv .* (contmod.p + contmod.dp)]

    @time begin
        for t in 1:Ndt
            #x = A \ (x + B)
            gmres!(x, A, x + B)
            if (mod(t, interval) == 0)
                thetas[:, Int64(t / interval)+1] = x[1:N]
                omegas[:, Int64(t / interval)+1] = x[N+1:end]
                ts[Int64(t / interval)+1] = t * dt
                println("NIter: $t Avg. Omega: $(sum(omegas[:, Int64(t/interval) + 1]) / N)")
            end
        end
    end
    return ts, thetas, omegas
end
