export get_params, update_model!, update_model_dm!, create_model

using SparseArrays


function create_model(
    mesh::ContGridMod.Mesh;
    filename::String = "",
    m_min::Real = 0.1,
    d_min::Real = 0.1,
    b_min::Real = 50,
)
    m = m_min * ones(mesh.Nedge)
    d = d_min * ones(mesh.Nedge)
    b = b_min * ones(mesh.Nedge)
    p = zeros(mesh.Nnode)
    th = zeros(mesh.Nnode)
    dp = zeros(mesh.Nnode)
    id_slack = 1
    
    if filename != ""
        dm = load_discrete_model(filename, mesh.scale_factor)
        id_slack, p, dp, b, m, d, th = get_params(mesh, dm,
            m_min = m_min, d_min = d_min, b_min = b_min)
    end
    
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


function get_params(
    mesh::ContGridMod.Mesh,
    dm::ContGridMod.DiscModel;
    m_min::Real = 0.1,
    d_min::Real = 0.1,
    b_min::Real = 50,
)    
    # load the discrete model
    
    m = zeros(mesh.Nnode)
    d = zeros(mesh.Nnode)
    pl = zeros(mesh.Nnode)
    pg = zeros(mesh.Nnode)
    th = zeros(mesh.Nnode)
    dp = zeros(mesh.Nnode)
    
    is_prod = dm.p_gen .> 0.0
    
    for (i, v) in enumerate(mesh.cell_vertices)
        is_in = in_polygon(dm.coord, v)
        d[i] = sum(dm.d_load[is_in])
        m[i] = sum(dm.m_gen[is_in[dm.id_gen] .* is_prod])
        d[i] = sum(dm.d_gen[is_in[dm.id_gen] .* is_prod])
        pg[i] = sum(dm.p_gen[is_in[dm.id_gen] .* is_prod])
        pl[i] = sum(dm.p_load[is_in])
        th[i] = sum(dm.th[is_in]) / sum(is_in)
    end
    
    p = pg - pl 
    p = p .- sum(p) / length(p)
    m = max.(m, m_min)
    d = max.(d, d_min)
    
    id_slack = find_nearest(dm.coord[dm.id_slack], mesh.node_coord)
    b = b_min * ones(mesh.Nedge) # TODO change that
    
    return id_slack, p, dp, b, m, d, th
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
