export get_params

using HDF5
using Plots
using SparseArrays

function get_params(
    mesh::Mesh,
    scale_factor::Float64,
    dataname::String;
    savename::String = "",
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
    dm = load_discrete_model(dataname, scale_factor)
    mesh_coord = mesh.coord[mesh.isgrid,:]
    idmesh = findall(mesh.isgrid)
    idin = findall(mesh.isinside)
    N = size(mesh.coord, 1)
    m = zeros(N)
    d = zeros(N)
    pl = zeros(N)
    pg = zeros(N)
    bx = zeros(N)
    by = zeros(N)
    @time begin
    Threads.@threads for g in 1:length(dm.idgen)
        if(dm.gen[g] > 0.0)
            k = argmin((mesh_coord[:, 1] .- dm.coord[dm.idgen[g],1]).^2 +
                (mesh_coord[:, 2] .- dm.coord[dm.idgen[g],2]).^2)
            m[idmesh[k]] += dm.mg[g]
            d[idmesh[k]] += dm.dg[g]
            pg[idmesh[k]] += dm.gen[g]
        end
    end
    
    Threads.@threads for l in 1:length(dm.dl)
        k = argmin((mesh_coord[:, 1] .- dm.coord[l,1]).^2 +
            (mesh_coord[:, 2] .- dm.coord[l,2]).^2)
        d[idmesh[k]] += dm.dl[l]
        pl[idmesh[k]] += dm.load[l]
    end

    bm = sum(dm.bline) / length(dm.bline)
    Threads.@threads for l in 1:size(dm.idb,1)
        x2 = dm.coord[dm.idb[l,2], 2]
        x1 = dm.coord[dm.idb[l,1], 2]
        y2 = dm.coord[dm.idb[l,2], 1]
        y1 = dm.coord[dm.idb[l,1], 1] 
        dx_l = x2 - x1
        dy_l = y2 - y1
        ds2 = (dy_l^2 + dx_l^2)
        phi = atan(dy_l, dx_l)
        if(dx_l != 0 && dy_l != 0) # if not a transformer
            y = mesh_coord[:, 1] # it's "less precise" but way simpler to implement 
            x = mesh_coord[:, 2] # using node locations instead of centers of lines, 
            beta = (dx_l .* (y1 .- y) + dy_l .* (x .- x1)) ./ ds2
            alpha = (dx_l .* (x .- x1) + dy_l .* (y .- y1)) ./ ds2
            in_seg = (0 .< alpha) .& (alpha .< 1)
            dist = abs.(beta) .* in_seg .* sqrt(ds2) + # if close to the segment
                .!in_seg .* min.(sqrt.((x .- x1).^2 + (y .- y1).^2), # if close to the ends
                sqrt.((x .- x2).^2 + (y .- y2).^2))
            if(dm.bline[l] < 2*bm) # if b is too large discard it
                bx[idmesh[dist .< dmax]] .+= dm.bline[l] * abs(cos(phi)) * dx^2 * patch
                by[idmesh[dist .< dmax]] .+= dm.bline[l] * abs(sin(phi)) * dx^2 * patch
            end
        end
    end
    end
    
    bx[mesh.isgrid] .= max.(bx[mesh.isgrid], bmin)
    by[mesh.isgrid] .= max.(by[mesh.isgrid], bmin)
    
    @time begin
        m = heat_diff(mesh, m, Niter = Niter, tau = tau)
        d = heat_diff(mesh, d, Niter = Niter, tau = tau)
        pl = heat_diff(mesh, pl, Niter = Niter, tau = tau)
        pg = heat_diff(mesh, pg, Niter = Niter, tau = tau)
        bx = heat_diff(mesh, bx, Niter = Niter, tau = tau)
        by = heat_diff(mesh, by, Niter = Niter, tau = tau)
    end

    # asign minimal values to the quantities
    m[mesh.isgrid] .= max.(m[mesh.isgrid], min_factor * sum(m) / sum(mesh.isgrid))
    d[mesh.isgrid] .= max.(d[mesh.isgrid], min_factor * sum(d) / sum(mesh.isgrid))
    #bx[isgrid] .= max.(bx[isgrid], min_factor * sum(bx) / sum(isgrid))
    #by[isgrid] .= max.(by[isgrid], min_factor * sum(by) / sum(isgrid))
    
    # ensure that the integrals of the different quantities over
    # the medium is equivalent to their sum over the discrete model
    m = (sum(dm.mg) / sum(m) / dx^2) .* m 
    d = ((sum(dm.dg) + sum(dm.dl)) / sum(d) / dx^2) .* d
    pl = (sum(dm.load) / sum(pl) / dx^2) .* pl
    pg = (sum(dm.gen) / sum(pg) / dx^2) .* pg
    
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
    
    for k in 1:mesh.Nx*mesh.Ny
        if(mesh.isinside[k])
            append!(id1, k)
            append!(id2, k)
            append!(v, - (bx[k] +
                bx[k-mesh.Ny] +
                by[k] +
                by[k-1])
                )
            append!(id1, k)
            append!(id2, k-mesh.Ny)
            append!(v, bx[k-mesh.Ny])
            append!(id1, k)
            append!(id2, k+mesh.Ny)
            append!(v, bx[k])
            append!(id1, k)
            append!(id2, k-1)
            append!(v, by[k-1])
            append!(id1, k)
            append!(id2, k+1)
            append!(v, by[k])
        end
    end
    
    for id in 1:size(mesh.n, 1)
        k = (Int64(mesh.n[id, 2]) - 1) * mesh.Ny + Int64(mesh.n[id, 1])
        ny = mesh.n[id, 3]
        nx = mesh.n[id, 4] 
        etamx = 1 - nx/2 - nx^2/2
        etapx = 1 + nx/2 - nx^2/2
        etamy = 1 - ny/2 - ny^2/2
        etapy = 1 + ny/2 - ny^2/2
        append!(id1, k)
        append!(id2, k)
        append!(v, - (etamx * bx[k] +
            etapx * bx[k-mesh.Ny] +
            etamy * by[k] +
            etapy * by[k-1]))
        append!(id1, k)
        append!(id2, k-mesh.Ny)
        append!(v, etapx * bx[k-mesh.Ny])
        append!(id1, k)
        append!(id2, k+mesh.Ny)
        append!(v, etamx * bx[k])
        append!(id1, k)
        append!(id2, k-1)
        append!(v, etapy * by[k-1])   
        append!(id1, k)
        append!(id2, k+1)
        append!(v, etamy * by[k])
    end
    
    xi = sparse(id1, id2, v, mesh.Ny * mesh.Nx, mesh.Ny * mesh.Nx)
    minv = m.^(-1)
    minv = minv[mesh.isgrid]
    gamma = d[mesh.isgrid] .* minv
    xi = SparseMatrixCSC{Float64, Int64}(xi[mesh.isgrid, mesh.isgrid])
    
    p = pg - pl 
    p = p[mesh.isgrid] .- sum(p[mesh.isgrid]) / sum(mesh.isgrid)
    
    cm = ContModel(
        mesh.Nx,
        mesh.Ny,
        mesh.coord,
        mesh.isinside,
        mesh.isborder,
        mesh.isgrid,
        mesh.yrange,
        mesh.xrange,
        mesh.n,
        mesh.dx,
        minv,
        gamma,
        p,
        xi,
        bx,
        by,
        m,
        d,
        zeros(size(d)))
    return cm
end


function heat_diff(
    mesh::mesh,
    q::Array{Float64,1};
    Niter::Int64 = 5000,
    tau::Float64 = 0.001,
)
    # parameters over the lattice with tau=kappa*dt

    id1 = Array{Int64,1}()
    id2 = Array{Int64,1}()
    v = Array{Float64,1}()
    for k in 1:length(q)
        if(mesh.isinside[k])
            append!(id1, k)
            append!(id2, k)
            append!(v, -4.0)
            append!(id1, k)
            append!(id2, k-mesh.Ny)
            append!(v, 1.0)
            append!(id1, k)
            append!(id2, k+mesh.Ny)
            append!(v, 1.0)
            append!(id1, k)
            append!(id2, k-1)
            append!(v, 1.0)
            append!(id1, k)
            append!(id2, k+1)
            append!(v, 1.0)            
        end
    end
    
    for id in 1:size(mesh.n, 1)
        k = (Int64(mesh.n[id, 2]) - 1) * mesh.Ny + Int64(mesh.n[id, 1])
        ny = mesh.n[id, 3]
        nx = mesh.n[id, 4] 
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
        append!(id2, k-mesh.Ny)
        append!(v, etapx)
        append!(id1, k)
        append!(id2, k+mesh.Ny)
        append!(v, etamx)
        append!(id1, k)
        append!(id2, k-1)
        append!(v, etapy)   
        append!(id1, k)
        append!(id2, k+1)
        append!(v, etamy)
    end
    
    N = mesh.Ny * mesh.Nx
    xi = sparse(id1, id2, v, N, N)
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
    coord = alberts_projection( data["bus_coord"] ./ (180 / pi) )
    coord = coord[:,[2,1]] / scaling_factor
    dm = DiscModel(
        vec(data["gen_inertia"]),
        vec(data["gen_prim_ctrl"]),
        Int64.(vec(data["gen"][:, 1])),
        coord,
        vec(data["load_freq_coef"]),
        Int64.(data["branch"][:, 1:2]),
        1.0 ./ data["branch"][:, 4],
        vec(data["bus"][:, 3]) / 100.0,
        vec(data["bus"][:, 9]) / 180.0 * pi,
        vec(data["gen"][:, 2]) / 100.0,
        vec(data["gen"][:, 9]) / 100.0,
        length(data["load_freq_coef"])
        )
    return dm
end


function update_params!(
    contmod::ContModel
)
    # update xi, minv and gamma

    id1 = Array{Int64,1}()
    id2 = Array{Int64,1}()
    v = Array{Float64,1}()
    
    for k in 1:contmod.Nx*contmod.Ny
        if(contmod.isinside[k])
            append!(id1, k)
            append!(id2, k)
            append!(v, - (contmod.bx[k] +
                contmod.bx[k-contmod.Ny] +
                contmod.by[k] +
                contmod.by[k-1])
                )
            append!(id1, k)
            append!(id2, k-contmod.Ny)
            append!(v, contmod.bx[k-contmod.Ny])
            append!(id1, k)
            append!(id2, k+contmod.Ny)
            append!(v, contmod.bx[k])
            append!(id1, k)
            append!(id2, k-1)
            append!(v, contmod.by[k-1])
            append!(id1, k)
            append!(id2, k+1)
            append!(v, contmod.by[k])
        end
    end
    
    for id in 1:size(contmod.n, 1)
        k = (Int64(contmod.n[id, 2]) - 1) * contmod.Ny + Int64(contmod.n[id, 1])
        ny = contmod.n[id, 3]
        nx = contmod.n[id, 4] 
        etamx = 1 - nx/2 - nx^2/2
        etapx = 1 + nx/2 - nx^2/2
        etamy = 1 - ny/2 - ny^2/2
        etapy = 1 + ny/2 - ny^2/2
        append!(id1, k)
        append!(id2, k)
        append!(v, - (etamx * contmod.bx[k] +
            etapx * contmod.bx[k-contmod.Ny] +
            etamy * contmod.by[k] +
            etapy * contmod.by[k-1]))
        append!(id1, k)
        append!(id2, k-contmod.Ny)
        append!(v, etapx * contmod.bx[k-contmod.Ny])
        append!(id1, k)
        append!(id2, k+contmod.Ny)
        append!(v, etamx * contmod.bx[k])
        append!(id1, k)
        append!(id2, k-1)
        append!(v, etapy * contmod.by[k-1])   
        append!(id1, k)
        append!(id2, k+1)
        append!(v, etamy * contmod.by[k])
    end
    
    xi = sparse(id1, id2, v, contmod.Ny * contmod.Nx, contmod.Ny * contmod.Nx)
    minv = contmod.m.^(-1)
    contmod.minv = minv[contmod.isgrid]
    contmod.gamma = contmod.d[contmod.isgrid] .* minv[contmod.isgrid]
    contmod.xi = SparseMatrixCSC{Float64, Int64}(xi[contmod.isgrid, contmod.isgrid])
end
