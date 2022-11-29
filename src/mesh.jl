export get_triangulated_mesh

using Triangulate


function get_triangulated_mesh(
    filename::String;
    max_area::Real = 0.5
)
    border, scale_factor = import_border(filename)

    # use Triangulate.jl to discretize the area
    seg = [collect(1:size(border,1)) [collect(2:size(border,1));1]]
    triin = Triangulate.TriangulateIO()
    triin.pointlist = Matrix{Cdouble}(border')
    triin.segmentlist = Matrix{Cint}(seg')
    triin.segmentmarkerlist = Vector{Int32}(collect(1:size(border,1)))
    (triout, vorout) = triangulate("pea$(max_area)DQ", triin)

    # define the center of mass of each triangle as a node
    point = triout.pointlist
    x, y = point[1,:], point[2,:]
    node = eachrow(triout.trianglelist') .|> id -> (sum(x[id])/3, sum(y[id])/3)
    area = eachrow(triout.trianglelist') .|> id -> 1/2*(x[id[1]]*(y[id[2]] - y[id[3]]) +
        x[id[2]]*(y[id[3]] - y[id[1]]) + x[id[3]]*(y[id[1]] - y[id[2]])) 
    cell = eachrow(triout.trianglelist') .|> id ->
        [(x[id[1]],y[id[1]]) , (x[id[2]],y[id[2]]), (x[id[3]],y[id[3]])]
    
    # find neighbours
    edge_list = Tuple{Int64,Int64}[]
    for i = 1:maximum(triout.trianglelist)
        list = findall(vec(any(triout.trianglelist .== i, dims=1)))
        for i=1:length(list)-1
            for j=i+1:length(list)
                # if they share an edge
                if length(unique([triout.trianglelist[:,list[i]]; triout.trianglelist[:,list[j]]])) == 4
                    push!(edge_list, (list[i], list[j]))
                end
            end
        end
    end

    #clear the list from duplicates
    edge = Tuple{Int64,Int64}[]
    for e = edge_list
        if all(edge .|> ne -> ne != e)
            push!(edge, e)
        end
    end
    border = eachrow(border) .|> b -> (b[1], b[2])
    
    return Mesh(node, edge, area, cell, border, scale_factor, length(edge), length(node))
end


function find_node(
    m::Mesh,
    coord::Vector{Tuple{Float64, Float64}},
)
    to_mat()
    id = Int64[]
    for c in coord
        # TODO
        push!(id, argmin((m.node_coord[:,1] .- coord[i,1]).^2 +
            (m.node_coord[:,2] .- coord[i,2]).^2))
    end
    return id
end
