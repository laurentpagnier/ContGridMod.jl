function add_local_disturbance!(
    contmod::ContModel,
    location::Array{Float64,1},
    dP::Float64,
    sigma::Float64,
    dx::Float64
)
    grid_coord = contmod.coord[contmod.isgrid,:]
    dp = zeros(size(grid_coord,1))
    for k = 1:size(grid_coord,1)
        dp[k] = exp(-((grid_coord[k,1] - location[1])^2 + (grid_coord[k,2] - location[2])^2)/2/sigma^2)
    end
    contmod.p += dP .* dp ./ sum(dp) / dx^2
    return dP .* dp ./ sum(dp) / dx^2
end



function add_db(
    isinside::BitVector,
    isgrid::BitVector,
    Ny::Int64,
    Nx::Int64,
    grid_coord::Array{Float64, 2},
    location::Array{Float64,1},
    dP::Float64,
    sigma::Float64,
    dx::Float64
)
    # flatten the m, d, p, bx and by matrices
    temp = local_disturbance(grid_coord, location, dP, sigma, dx)
    db = zeros(size(isinside))
    db[isgrid] = temp
    id1 = Array{Int64,1}()
    id2 = Array{Int64,1}()
    v = Array{Float64,1}()
    for k in 1:Nx*Ny
        if(isinside[k])
            append!(id1, k)
            append!(id2, k)
            append!(v, - (db[k] +
                db[k-Ny] +
                db[k] +
                db[k-1])
                )
            append!(id1, k)
            append!(id2, k-Ny)
            append!(v, db[k-Ny])
            append!(id1, k)
            append!(id2, k+Ny)
            append!(v, db[k])
            append!(id1, k)
            append!(id2, k-1)
            append!(v, db[k-1])
            append!(id1, k)
            append!(id2, k+1)
            append!(v, db[k])            
        end
    end

    for id in 1:size(n, 1)
        k = (Int64(n[id, 2]) - 1) * Ny + Int64(n[id, 1])
        ny = n[id, 3]
        nx = n[id, 4] 
        etamx = 1 - nx/2-nx^2/2
        etapx = 1 + nx/2-nx^2/2
        etamy = 1 - ny/2-ny^2/2
        etapy = 1 + ny/2-ny^2/2
        append!(id1, k)
        append!(id2, k)
        append!(v, - (etamx * db[k] +
            etapx * db[k-Ny] +
            etamy * db[k] +
            etapy * db[k-1]))
        append!(id1, k)
        append!(id2, k-Ny)
        append!(v, etapx * db[k-Ny])
        append!(id1, k)
        append!(id2, k+Ny)
        append!(v, etamx * db[k])
        append!(id1, k)
        append!(id2, k-1)
        append!(v, etapy * db[k-1])   
        append!(id1, k)
        append!(id2, k+1)
        append!(v, etamy * db[k])
    end
    
    xi = sparse(id1, id2, v, Ny * Nx, Ny * Nx)
    xi = SparseMatrixCSC{Float64, Int64}(xi[isgrid,isgrid])
    return xi, db
end
