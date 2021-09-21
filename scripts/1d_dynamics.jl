function perform_dyn_sim(
    b::Array{Float64, 1},
    p::Array{Float64, 1},
    m::Array{Float64, 1},
    d::Array{Float64, 1},
    th0::Array{Float64, 1};
    interval::Int64 = 100,
    Ndt::Int64 = 15000,
    dt = 0.0001,
    boundary::String = "third"
)
    
    minvflat = 1 ./ m
    gammaflat = d ./ m

    th = copy(th0)
    th_new = zeros(size(th0))
    th_old = copy(th0)

    ts = zeros(1 + Int64(ceil(Ndt/interval)))
    
    omegas = zeros(N, 1 + Int64(ceil(Ndt/interval)))
    thetas = zeros(N, 1 + Int64(ceil(Ndt/interval)))
    th = copy(th0)
    th_new = zeros(N)
    th_old = copy(th0)
    omegas[:,1] = (th - th_old) / dt;
    thetas[:,1] = th;
    println(sum(p)*dx)
    ts = zeros(1 + Int64(ceil(Ndt/interval)))

    chi = 1.0 .+ gammaflat*dt / 2.0
    xi = sparse([], [], Float64.([]), N, N)
    
    #Threads.@threads for k in 1:Nx*Ny
    for k in 2:N-1
        xi[k, k] = - b[k] - b[k-1]
        xi[k, k-1] = b[k-1]
        xi[k, k+1] = b[k]
    end

    if(boundary == "fourth")
        # fourth order boundary condition
        println("fourth order")
        xi[1,1] = - 3.5 * b[1]
        xi[1,2] = 7.0 * b[1]
        xi[1,3] = -0.5 * b[1]
        xi[N,N] = - 3.5 * b[N-1]
        xi[N,N-1] = 4.0 * b[N-1]
        xi[N,N-2] = -0.5 * b[N-1]

    elseif(boundary == "fifth")
        println("fifth order")
        xi[1,1] = - 16.0 / 3.0 * b[1]
        xi[1,2] =  7.0 * b[1]
        xi[1,3] = -2.0 * b[1]
        xi[1,4] = 1 / 3.0 * b[1]
        xi[N,N] = - 16.0 / 3.0 * b[N-1]
        xi[N,N-1] = 7.0 * b[N-1]
        xi[N,N-2] = -2.0 * b[N-1]
        xi[N,N-3] = 1 / 3.0 * b[N-1]
    else
        # thrid order boundary condition
        println("third order")
        xi[1,1] = - 2 * b[1]
        xi[1,2] = 2 * b[1]
        xi[N,N] = - 2 * b[N-1]
        xi[N,N-1] = 2 * b[N-1]
    end
    xi = SparseMatrixCSC{Float64, Int64}(xi)

    chi = 1 ./ (1 .+ gammaflat*dt/2)

    A = (2 * chi - dt^2/dx^2 * minvflat)
    B = (1 .- gammaflat * dt / 2) .* chi
    C = dt^2 * chi .* minvflat .* p

    @time begin
        for t in 1:Ndt
            th_new = 2 * chi .* th +
                dt^2 / dx^2 .* minvflat .* chi .* (xi * th) -
                B .* th_old + C
        
            if(mod(t, interval) == 0)
                println("NIter: ", t)
                omegas[:,Int64(t/interval) + 1] = (th_new - th) / dt
                thetas[:,Int64(t/interval) + 1] = th_new
                ts[Int64(t/interval) + 1] = t * dt
            end
            
            th_old = copy(th)
            th = copy(th_new)
        end
    end
    return ts, thetas, omegas
end
