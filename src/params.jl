export get_params, update_model!, update_model_dm!

using SparseArrays

function get_params(
    mesh::ContGridMod.Mesh,
    dataname::String;
    m_min::Real = 0.1,
    d_min::Real = 0.1,
    b_min::Real = 50,
)    
    # load the discrete model
    dm = ContGridMod.load_discrete_model(dataname, mesh.scale_factor)
    m = zeros(mesh.Nnode)
    d = zeros(mesh.Nnode)
    pl = zeros(mesh.Nnode)
    pg = zeros(mesh.Nnode)
    th = zeros(mesh.Nnode)
    dp = zeros(mesh.Nnode)
    
    is_prod = dm.p_gen .> 0.0
    
    for (i, v) in enumerate(mesh.cell_vertices)
        is_in = ContGridMod.in_polygon(dm.coord, v)
        d[i] += sum(dm.d_load[is_in])
        m[i] += sum(dm.m_gen[is_in[dm.id_gen] .* is_prod])
        d[i] += sum(dm.d_gen[is_in[dm.id_gen] .* is_prod])
        pg[i] += sum(dm.p_gen[is_in[dm.id_gen] .* is_prod])
        pl[i] += sum(dm.p_load[is_in])
    end
    
    p = pg - pl 
    p = p .- sum(p) / length(p)
    m = max.(m, m_min)
    d = max.(d, d_min)
    
    id_slack = find_nearest(dm.coord[dm.id_slack], mesh.node_coord)
    b = 200 * ones(mesh.Nedge)
    
    return ContGridMod.ContModel(
        mesh,
        id_slack,
        p,
        dp,
        b,
        m,
        d,
        th,
    )
end


function heat_diff(
    mesh::Mesh,
    q::Vector{Float64};
    Niter::Int64 = 5000,
    tau::Float64 = 1.0,
)
    # parameters over the lattice with tau=kappa*dt
    # incidence matrix
    temp = to_mat(mesh.edge_list)
    inc = sparse([temp[:,1]; temp[:,2]],
        [1:mesh.Nedge; 1:mesh.Nedge], [-ones(mesh.Nedge); ones(mesh.Nedge)])
    # - Laplacian
    L = inc * inc' 
    
    # Crank-Nicolson
    I = sparse(1:mesh.Nnode, 1:mesh.Nnode, ones(mesh.Nnode))
    A = 2.0 .* I + tau .* L
    B = 2.0 .* I - tau .* L
    for _ in 1:Niter
        q = A \ (B * q)
    end

    return q
end
