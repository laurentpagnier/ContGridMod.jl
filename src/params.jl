export get_params, update_params!, load_discrete_model, update_model!, update_susceptance!

using HDF5
using Plots
using SparseArrays

function get_params(
    mesh::Mesh,
    scale_factor::Float64,
    dataname::String;
    dmax::Float64 = 100.0,
    Niter::Int64 = 5000,
    tau::Float64 = 0.001,
    patch::Float64 = 0.001,
    min_factor::Float64 = 0.1,
)
    # this function uses the heat equation to "diffuse" the discrete
    # parameters
    
    # load the discrete model
    dm = load_discrete_model(dataname, scale_factor)
    #idin = findall(mesh.isinside)
    N = size(mesh.coord, 1)
    m= zeros(N)
    d = zeros(N)
    pl = zeros(N)
    pg = zeros(N)
    b = zeros(N)
    @time begin
    Threads.@threads for g in 1:length(dm.idgen)
        if(dm.gen[g] > 0.0)
            k = argmin((mesh.coord[:, 1] .- dm.coord[dm.idgen[g],1]).^2 +
                (mesh.coord[:, 2] .- dm.coord[dm.idgen[g],2]).^2)
            m[k] += dm.mg[g]
            d[k] += dm.dg[g]
            pg[k] += dm.gen[g]
        end
    end
    
    Threads.@threads for l in 1:length(dm.dl)
        k = argmin((mesh.coord[:, 1] .- dm.coord[l,1]).^2 +
            (mesh.coord[:, 2] .- dm.coord[l,2]).^2)
        d[k] += dm.dl[l]
        pl[k] += dm.load[l]
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
        if(dx_l != 0 && dy_l != 0) # if not a transformer
            y = mesh.coord[:, 1] # it's "less precise" but way simpler to implement 
            x = mesh.coord[:, 2] # using node locations instead of centers of lines, 
            beta = (dx_l .* (y1 .- y) + dy_l .* (x .- x1)) ./ ds2
            alpha = (dx_l .* (x .- x1) + dy_l .* (y .- y1)) ./ ds2
            in_seg = (0 .< alpha) .& (alpha .< 1)
            dist = abs.(beta) .* in_seg .* sqrt(ds2) + # if close to the segment
                .!in_seg .* min.(sqrt.((x .- x1).^2 + (y .- y1).^2), # if close to the ends
                sqrt.((x .- x2).^2 + (y .- y2).^2))
            if(dm.bline[l] < 2*bm) # if b is too large discard it
                b[dist .< dmax] .+= dm.bline[l] * mesh.dx^2 * patch
            end
        end
    end
    
    
    end # time
    
    n = size(mesh.inc_mat, 1)
    N = size(mesh.coord, 1)
    
    @time begin
        m = heat_diff(mesh, m, Niter = Niter, tau = tau)
        d = heat_diff(mesh, d, Niter = Niter, tau = tau)
        pl = heat_diff(mesh, pl, Niter = Niter, tau = tau)
        pg = heat_diff(mesh, pg, Niter = Niter, tau = tau)
        b = heat_diff(mesh, b, Niter = Niter, tau = tau)
    end

    # asign minimal values to the quantities
    m .= max.(m, min_factor * sum(m) / length(m))
    d .= max.(d, min_factor * sum(d) / length(m))
    
    # ensure that the integrals of the different quantities over
    # the medium is equivalent to their sum over the discrete model
    m = (sum(dm.mg) / sum(m) /mesh.dx^2) .* m 
    d = ((sum(dm.dg) + sum(dm.dl)) / sum(d) / mesh.dx^2) .* d
    pl = (sum(dm.load) / sum(pl) / mesh.dx^2) .* pl
    pg = (sum(dm.gen) / sum(pg) / mesh.dx^2) .* pg

    # ensure that the integral of the power injection is 0, i.e
    # that generation matches the demand.

    temp = sparse([mesh.inc_mat[:,1]; mesh.inc_mat[:,2]], [1:n; 1:n], 0.5*[ones(n); ones(n)])
    b = temp' * b
    temp = sparse([mesh.inc_mat[:,1]; mesh.inc_mat[:,2]], [1:n; 1:n], [-ones(n); ones(n)])
    xi = - temp * (b.* temp')
    minv = m.^(-1)
    gamma = d.* minv
    
    p = pg - pl 
    p = p .- sum(p) / length(p)
    
    return cm = ContModel(
        mesh,
        minv,
        gamma,
        p,
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


function heat_diff(
    mesh::Mesh,
    q::Array{Float64,1};
    Niter::Int64 = 5000,
    tau::Float64 = 0.001,
)
    # parameters over the lattice with tau=kappa*dt
    n = size(mesh.inc_mat, 1)
    N = size(mesh.coord, 1)
    # incidence matrix
    temp = sparse([mesh.inc_mat[:,1]; mesh.inc_mat[:,2]], [1:n; 1:n], [-ones(n); ones(n)])
    # - Laplacian
    xi = -temp * temp' 
    
    # Crank-Nicolson
    I = sparse(1:N, 1:N, ones(N))
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
    coord = albers_projection( data["bus_coord"] ./ (180 / pi) )
    coord ./= scaling_factor
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
    minv = contmod.m.^(-1)
    gamma = contmod.d.* minv
    contmod.minv = minv
    contmod.gamma = contmod.d .* minv
    n = size(contmod.mesh.inc_mat, 1)
    temp = sparse([contmod.mesh.inc_mat[:,1]; contmod.mesh.inc_mat[:,2]], [1:n; 1:n], [-ones(n); ones(n)])
    contmod.xi = - temp * (contmod.b .* temp')
end

function update_model!(
    contmod::ContModel,
    dm::DiscModel
)
    # Update the dynamical parameters of the model from a discrete model
    N = size(contmod.mesh.coord, 1)

    m = zeros(N)
    d = zeros(N)
    pl = zeros(N)
    pg = zeros(N)

    Threads.@threads for g in 1:length(dm.idgen)
        if(dm.gen[g] > 0.0)
            k = argmin((contmod.mesh.coord[:, 1] .- dm.coord[dm.idgen[g],1]).^2 +
                (contmod.mesh.coord[:, 2] .- dm.coord[dm.idgen[g],2]).^2)
            m[k] += dm.mg[g]
            d[k] += dm.dg[g]
            pg[k] += dm.gen[g]
        end
    end

    Threads.@threads for l in 1:length(dm.dl)
        k = argmin((contmod.mesh.coord[:, 1] .- dm.coord[l,1]).^2 +
            (contmod.mesh.coord[:, 2] .- dm.coord[l,2]).^2)
        d[k] += dm.dl[l]
        pl[k] += dm.load[l]
    end
    @time begin
        m = heat_diff(contmod.mesh, m, Niter = contmod.Niter, tau = contmod.tau)
        d = heat_diff(contmod.mesh, d, Niter = contmod.Niter, tau = contmod.tau)
        pl = heat_diff(contmod.mesh, pl, Niter = contmod.Niter, tau = contmod.tau)
        pg = heat_diff(contmod.mesh, pg, Niter = contmod.Niter, tau = contmod.tau)
    end

    # asign minimal values to the quantities
    m .= max.(m, contmod.min_factor * sum(m) / length(m))
    d .= max.(d, contmod.min_factor * sum(d) / length(m))
    
    # ensure that the integrals of the different quantities over
    # the medium is equivalent to their sum over the discrete model
    m = (sum(dm.mg) / sum(m) / contmod.mesh.dx^2) .* m 
    d = ((sum(dm.dg) + sum(dm.dl)) / sum(d) / contmod.mesh.dx^2) .* d
    pl = (sum(dm.load) / sum(pl) / contmod.mesh.dx^2) .* pl
    pg = (sum(dm.gen) / sum(pg) / contmod.mesh.dx^2) .* pg

    # ensure that the integral of the power injection is 0, i.e
    # that generation matches the demand.
    p = pg - pl
    p = p .- sum(p) / length(p) 

    minv = m.^(-1)    
    # Update the paramters
    contmod.minv = minv
    contmod.gamma = d .* minv
    contmod.m = m
    contmod.d = d
    contmod.p = p
end

function update_susceptance!(cm::ContModel, b::Vector{Float64})
    n = size(cm.mesh.inc_mat, 1)
    incMat = sparse([cm.mesh.inc_mat[:,1]; cm.mesh.inc_mat[:,2]], [1:n; 1:n], [-ones(n); ones(n)])
    cm.b = b
    cm.xi = -incMat * (b .* incMat')
end
