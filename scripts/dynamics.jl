using SparseArrays
using LinearAlgebra
using IterativeSolvers

function perform_dyn_sim(
    contmod::ContModel;
    interval::Int64,
    Ndt::Int64,
    dt::Float64,
    method::String="crank-nicolson"
)

    if(method == "crank-nicolson")
        return ts, thetas, omegas = perform_dyn_sim_crank_nicolson(
            contmod; interval, Ndt, dt)  
    elseif(method == "backward-euler")
        return ts, thetas, omegas = perform_dyn_sim_backward_euler(
            contmod; interval, Ndt, dt)    
    elseif(method == "forward")
        return ts, thetas, omegas = perform_dyn_sim_forward(
            contmod; interval, Ndt, dt)
    else
        println("Method not found, '", method, "' is not implemented.")
    end

end


function perform_dyn_sim_forward(
    contmod::ContModel;
    interval::Int64 = 100,
    Ndt::Int64 = 15000,
    dt::Float64 = 0.0001
)
    println("Total time: ", dt * Ndt)
    N = sum(contmod.isgrid)
    M = 1 + Int64(ceil(Ndt/interval))
    omegas = zeros(N, M)
    thetas = zeros(N, M)
    th_new = zeros(N)
    ts = zeros(M)
    
    th_old = copy(contmod.th[contmod.isgrid])
    th = copy(contmod.th[contmod.isgrid])
    omegas[:,1] = zeros(size(th))
    thetas[:,1] = copy(contmod.th[contmod.isgrid])  

    chi = 1 ./ (1 .+ contmod.gamma*dt/2)

    A = (2 * chi - dt^2/dx^2 * contmod.minv)
    B = (1 .- contmod.gamma * dt / 2) .* chi
    C = dt^2 * chi .* contmod.minv .* contmod.p
    
    @time begin
        for t in 1:Ndt
            th_new = 2 * chi .* th +
                dt^2 / dx^2 .* contmod.minv .* chi .* (xi * th) -
                B .* th_old + C
            if(mod(t,interval) == 0)
                omegas[:,Int64(t/interval) + 1] = (th_new-th) / dt
                thetas[:,Int64(t/interval) + 1] = th_new
                ts[Int64(t/interval) + 1] = t * dt
                println("NIter: ", t, " Avg. Omega: ", sum(omegas[:, Int64(t/interval) + 1])/sum(contmod.isgrid))
            end
            th_old = copy(th)
            th = copy(th_new)
        end
    end

    return ts, thetas, omegas
end


function perform_dyn_sim_crank_nicolson(
    contmod::ContModel;
    interval::Int64 = 10,
    Ndt::Int64 = 1000,
    dt::Float64 = 0.05
)
    println("Total time: ", dt * Ndt)
    N = sum(contmod.isgrid)
    M = 1 + Int64(ceil(Ndt/interval))
    omegas = zeros(N, M)
    thetas = zeros(N, M)
    ts = zeros(M)
    
    x = [contmod.th[contmod.isgrid]; zeros(N)]
    omegas[:, 1] = zeros(N)
    thetas[:, 1] = copy(contmod.th[contmod.isgrid])  

    ts = zeros(1 + Int64(ceil(Ndt/interval)))

    I = sparse(1:N, 1:N, ones(N))
    A = [I -dt / 2 * I;
        - dt / 2 / dx^2 * sparse(1:N, 1:N, contmod.minv) * contmod.xi (I + dt/2 * sparse(1:N, 1:N, contmod.gamma))]
    B = [I dt / 2 * I;
         dt / 2 / dx^2 * sparse(1:N, 1:N, contmod.minv) * contmod.xi (I - dt/2 * sparse(1:N, 1:N, contmod.gamma))]
    C = [zeros(N); dt * sparse(1:N, 1:N, contmod.minv) * contmod.p]

    @time begin
        for t in 1:Ndt
            #x = A \ (B * x + C) # way slower when dx -> 0
            gmres!(x, A , B * x + C)
            if(mod(t,interval) == 0)
                thetas[:,Int64(t/interval) + 1] = x[1:N]
                omegas[:,Int64(t/interval) + 1] = x[N+1:end]
                ts[Int64(t/interval) + 1] = t*dt
                println("NIter: ", t, " Avg. Omega: ", sum(omegas[:, Int64(t/interval) + 1]) / sum(contmod.isgrid))
            end
        end
    end
    return ts, thetas, omegas
end


function perform_dyn_sim_backward_euler(
    contmod::ContModel;
    interval::Int64 = 10,
    Ndt::Int64 = 1000,
    dt::Float64 = 0.05
)
    println("Total time: ", dt * Ndt)
    N = sum(contmod.isgrid)
    M = 1 + Int64(ceil(Ndt/interval))
    omegas = zeros(N, M)
    thetas = zeros(N, M)
    ts = zeros(M)
    
    x = [contmod.th[contmod.isgrid]; zeros(N)]
    omegas[:,1] = zeros(N)
    thetas[:,1] = copy(contmod.th[contmod.isgrid])  

    ts = zeros(1 + Int64(ceil(Ndt/interval)))

    I = sparse(1:N, 1:N, ones(N))
    A = [I -dt * I;
        - dt / dx^2 * sparse(1:N, 1:N, contmod.minv) * contmod.xi (I + dt * sparse(1:N, 1:N, contmod.gamma))]
    B = [zeros(N); dt * sparse(1:N, 1:N, contmod.minv) * contmod.p]

    @time begin
        for t in 1:Ndt
            #x = A \ (x + B)
            gmres!(x, A , x + B)
            if(mod(t,interval) == 0)
                thetas[:,Int64(t/interval) + 1] = x[1:N]
                omegas[:,Int64(t/interval) + 1] = x[N+1:end]
                ts[Int64(t/interval) + 1] = t*dt
                if(mod(t, 20) == 0)
                    println("NIter: ", t, " Avg. Omega: ", sum(omegas[:, Int64(t/interval) + 1])/sum(contmod.isgrid))
                end
            end
        end
    end
    return ts, thetas, omegas
end
