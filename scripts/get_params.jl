using HDF5
using Plots

function get_params(
    isinside::BitVector,
    isgrid::BitVector,
    n::Array{Float64, 2},
    Ny::Int64,
    Nx::Int64,
    dx::Float64,
    lat_coord::Array{Float64, 2},
    scale_factor::Float64,
    dataname::String,
    savename::String;
    def_val::Float64 = 1E-2,
    dmax::Float64 = 100.0,
    Niter::Int64 = 5000,
    tau::Float64 = 0.001,
    patch::Float64 = 0.001,
    min_factor::Float64 = 0.1,
    bmin::Float64 = 1.0
)
    # this function uses the heat equation to "diffuse" the discrete
    # parameters
    
    # load the discrete model
    gen, dem, bline, idb, idgen, coord, mg, dg, dl, th = load_discrete_model(dataname, scale_factor)
    
    grid_coord = lat_coord[isgrid,:]
    idgrid = findall(isgrid)
    idin = findall(isinside)
    N = size(lat_coord, 1)
    m = zeros(N)
    d = zeros(N)
    pl = zeros(N)
    pg = zeros(N)
    bx = zeros(N)
    by = zeros(N)
    @time begin
    Threads.@threads for g in 1:length(idgen)
        if(gen[g] > 0.0)
            k = argmin((grid_coord[:, 1] .- coord[idgen[g],1]).^2 +
                (grid_coord[:, 2] .- coord[idgen[g],2]).^2)
            m[idgrid[k]] += mg[g]
            d[idgrid[k]] += dg[g]
            pg[idgrid[k]] += gen[g]
        end
    end
    
    Threads.@threads for l in 1:length(dl)
        k = argmin((grid_coord[:, 1] .- coord[l,1]).^2 +
            (grid_coord[:, 2] .- coord[l,2]).^2)
        d[idgrid[k]] += dl[l]
        pl[idgrid[k]] += dem[l]
    end

    bm = sum(bline) / length(bline)
    Threads.@threads for l in 1:size(idb,1)
        x2 = coord[idb[l,2], 2]
        x1 = coord[idb[l,1], 2]
        y2 = coord[idb[l,2], 1]
        y1 = coord[idb[l,1], 1] 
        dx_l = x2 - x1
        dy_l = y2 - y1
        ds2 = (dy_l^2 + dx_l^2)
        phi = atan(dy_l, dx_l)
        if(dx_l != 0 && dy_l != 0) # if not a transformer
            y = grid_coord[:, 1] # it's "less precise" but way simpler to implement 
            x = grid_coord[:, 2] # using node locations instead of centers of lines, 
            beta = (dx_l .* (y1 .- y) + dy_l .* (x .- x1)) ./ ds2
            alpha = (dx_l .* (x .- x1) + dy_l .* (y .- y1)) ./ ds2
            in_seg = (0 .< alpha) .& (alpha .< 1)
            dist = abs.(beta) .* in_seg .* sqrt(ds2) + # if close to the segment
                .!in_seg .* min.(sqrt.((x .- x1).^2 + (y .- y1).^2), # if close to the ends
                sqrt.((x .- x2).^2 + (y .- y2).^2))
            if(bline[l] < 2*bm)
                bx[idgrid[dist .< dmax]] .+= bline[l] * abs(cos(phi)) * dx^2 * patch
                by[idgrid[dist .< dmax]] .+= bline[l] * abs(sin(phi)) * dx^2 * patch
            end
        end
    end
    end
    
    bx[isgrid] .= max.(bx[isgrid], bmin)
    by[isgrid] .= max.(by[isgrid], bmin)
    
    @time begin
        m = heat_diff(isinside, n, m, dx, Ny, Nx, Niter = Niter, tau = tau)
        d = heat_diff(isinside, n, d, dx, Ny, Nx, Niter = Niter, tau = tau)
        pl = heat_diff(isinside, n, pl, dx, Ny, Nx, Niter = Niter, tau = tau)
        pg = heat_diff(isinside, n, pg, dx, Ny, Nx, Niter = Niter, tau = tau)
        bx = heat_diff(isinside, n, bx, dx, Ny, Nx, Niter = Niter, tau = tau)
        by = heat_diff(isinside, n, by, dx, Ny, Nx, Niter = Niter, tau = tau)
    end

    # asign minimal values to the quantities
    m[isgrid] .= max.(m[isgrid], min_factor * sum(m) / sum(isgrid))
    d[isgrid] .= max.(d[isgrid], min_factor * sum(d) / sum(isgrid))
    #bx[isgrid] .= max.(bx[isgrid], min_factor * sum(bx) / sum(isgrid))
    #by[isgrid] .= max.(by[isgrid], min_factor * sum(by) / sum(isgrid))
    
    # ensure that the integrals of the different quantities over
    # the medium is equivalent to their sum over the discrete model
    m = (sum(mg) / sum(m) / dx^2) .* m 
    d = ((sum(dg) + sum(dl)) / sum(d) / dx^2) .* d
    pl = (sum(dem) / sum(pl) / dx^2) .* pl
    pg = (sum(gen) / sum(pg) / dx^2) .* pg
    
    # save the quantities
    #fid = h5open(savename, "w")
    #write(fid, "bx", bx)
    #write(fid, "by", by)
    #write(fid, "m", m)
    #write(fid, "d", d)
    #write(fid, "pl", pl)
    #write(fid, "pg", pg)
    #write(fid, "xrange", xrange)
    #write(fid, "yrange", yrange)
    #close(fid)

    # ensure that the integral of the power injection is 0, i.e
    # that generation matches the demand.

    id1 = Array{Int64,1}()
    id2 = Array{Int64,1}()
    v = Array{Float64,1}()
    for k in 1:Nx*Ny
        if(isinside[k])
            append!(id1, k)
            append!(id2, k)
            append!(v, - (bx[k] +
                bx[k-Ny] +
                by[k] +
                by[k-1])
                )
            append!(id1, k)
            append!(id2, k-Ny)
            append!(v, bx[k-Ny])
            append!(id1, k)
            append!(id2, k+Ny)
            append!(v, bx[k])
            append!(id1, k)
            append!(id2, k-1)
            append!(v, by[k-1])
            append!(id1, k)
            append!(id2, k+1)
            append!(v, by[k])
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
        append!(v, - (etamx * bx[k] +
            etapx * bx[k-Ny] +
            etamy * by[k] +
            etapy * by[k-1]))
        append!(id1, k)
        append!(id2, k-Ny)
        append!(v, etapx * bx[k-Ny])
        append!(id1, k)
        append!(id2, k+Ny)
        append!(v, etamx * bx[k])
        append!(id1, k)
        append!(id2, k-1)
        append!(v, etapy * by[k-1])   
        append!(id1, k)
        append!(id2, k+1)
        append!(v, etamy * by[k])
    end
    
    xi = sparse(id1, id2, v, Ny * Nx, Ny * Nx)
    minv = m.^(-1)
    minv = minv[isgrid]
    gamma = d[isgrid] .* minv
    xi = SparseMatrixCSC{Float64, Int64}(xi[isgrid, isgrid])
    
    p = pg - pl 
    p = p[isgrid] .- sum(p[isgrid]) / sum(isgrid)
    return minv, gamma, p, xi, bx, by, m, d
end





function heat_diff(
    ininside::BitVector,
    n::Array{Float64,2},
    q::Array{Float64,1},
    dx::Float64,
    Ny::Int64,
    Nx::Int64;
    Niter::Int64 = 5000,
    tau::Float64 = 0.001,
)
    # parameters over the lattice with tau=kappa*dt

    id1 = Array{Int64,1}()
    id2 = Array{Int64,1}()
    v = Array{Float64,1}()
    for k in 1:length(q)
        if(ininside[k])
            append!(id1, k)
            append!(id2, k)
            append!(v, -4.0)
            append!(id1, k)
            append!(id2, k-Ny)
            append!(v, 1.0)
            append!(id1, k)
            append!(id2, k+Ny)
            append!(v, 1.0)
            append!(id1, k)
            append!(id2, k-1)
            append!(v, 1.0)
            append!(id1, k)
            append!(id2, k+1)
            append!(v, 1.0)            
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
        append!(v, - (etamx +
            etapx +
            etamy +
            etapy))
        append!(id1, k)
        append!(id2, k-Ny)
        append!(v, etapx)
        append!(id1, k)
        append!(id2, k+Ny)
        append!(v, etamx)
        append!(id1, k)
        append!(id2, k-1)
        append!(v, etapy)   
        append!(id1, k)
        append!(id2, k+1)
        append!(v, etamy)
    end

    xi = sparse(id1, id2, v, Ny * Nx, Ny * Nx)
    N = Ny * Nx
    I = sparse(1:N, 1:N, ones(N))
    A = 2.0 .* I - tau .* xi / dx^2
    B = 2.0 .* I + tau .* xi / dx^2
    for t in 1:Niter
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
    mg = vec(data["gen_inertia"])
    dg = vec(data["gen_prim_ctrl"])
    idgen = Int64.(vec(data["gen"][:, 1]))
    coord = alberts_projection( data["bus_coord"] ./ (180 / pi) )
    dl = vec(data["load_freq_coef"])
    idb = Int64.(data["branch"][:, 1:2])
    bline = 1 ./ data["branch"][:, 4]
    dem = vec(data["bus"][:, 3]) / 100.0
    th = vec(data["bus"][:, 9]) / 180.0 * pi
    gen = vec(data["gen"][:, 2]) / 100.0
    return gen, dem, bline, idb, idgen, coord[:,[2,1]] / scaling_factor, mg, dg, dl, th
end
