export get_params, update_model!, update_model_dm!, create_model,
    load_model, save_model

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
    dp(t) = zeros(mesh.Nnode)
    id_slack = 1
    
    if filename != ""
        dm = load_discrete_model(filename, mesh.scale_factor)
        id_slack, p, b, m, d, th = get_params(mesh, dm,
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


function save_model(
    filename::String,
    contmod::ContModel,
)
    fid = h5open(filename,"w")
    for f in fieldnames(ContGridMod.ContModel)
        v = getproperty(contmod, f)
        if typeof(v) == ContGridMod.Mesh
            fid["mesh/node_coord"] = reduce(vcat, v.node_coord .|> v -> [v[1] v[2]])
            fid["mesh/edge_list"] = reduce(vcat, v.edge_list .|> v -> [v[1] v[2]])
            fid["mesh/cell_area"] = v.cell_area
            for (i, w) in enumerate(v.cell_vertices)
                fid["mesh/cell_vertices/$i"] = reduce(vcat, w .|> v -> [v[1] v[2]])
            end
            fid["mesh/border"] = reduce(vcat, v.border .|> v -> [v[1] v[2]])
            fid["mesh/scale_factor"] = v.scale_factor
            fid["mesh/Nedge"] = v.Nedge
            fid["mesh/Nnode"] = v.Nnode
        elseif !(v isa Function)
            fid["$f"] = v
        end
    end
    close(fid)
end


function load_model(
    filename::String
)
    d = h5read(filename, "/")
    node_coord = eachrow(d["mesh"]["node_coord"]) .|> v -> (v[1], v[2])
    edge_list = eachrow(d["mesh"]["edge_list"]) .|> v -> (v[1], v[2])
    cell_area = d["mesh"]["cell_area"]
    cell_vertices = Vector{Tuple{Float64,Float64}}[]
    for i = 1:length(d["mesh"]["cell_vertices"])
        cv = eachrow(d["mesh"]["cell_vertices"]["$i"]) .|> v -> (v[1], v[2])
        push!(cell_vertices, cv)
    end
    border = d["mesh"]["border"]
    scale_factor = d["mesh"]["scale_factor"]
    Nedge = d["mesh"]["Nedge"]
    Nnode = d["mesh"]["Nnode"]

    border = eachrow(d["mesh"]["border"]) .|> v -> (v[1], v[2])
    mesh = ContGridMod.Mesh(node_coord, edge_list, cell_area,
        cell_vertices, border, scale_factor, Nedge, Nnode)

    id_slack = d["id_slack"]
    p = d["p"]
    b = d["b"]
    m = d["m"]
    th = d["th"]
    d = d["d"]

    dp(t) = zeros(Nnode)
    ContGridMod.ContModel(mesh, id_slack, p,
       dp, b, m, d, th)
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
    
    return id_slack, p, b, m, d, th
end
