module ContGridMod

using JSON
using SparseArrays

mutable struct DiscModel
    m_gen::Vector{Float64} # generator inertia
    d_gen::Vector{Float64} # generator damping
    id_gen::Vector{Int64} # generator inertia
    id_slack::Int64 # index of the slack bus
    coord::Vector{Tuple{Float64, Float64}} # coordinates in 
    d_load::Vector{Float64}
    line_list::Vector{Tuple{Int64,Int64}} # list of indices
    b::Vector{Float64} # susceptance
    p_load::Vector{Float64}
    th::Vector{Float64}
    p_gen::Vector{Float64}
    max_gen::Vector{Float64}
    Nbus::Int64
    Ngen::Int64
    Nline::Int64
end


mutable struct FreqProfile
    coord::Vector{Tuple{Float64, Float64}}
    time::Vector{Float64}
    freq::Vector{Vector{Float64}}
end


mutable struct Mesh
    node_coord::Vector{Tuple{Float64,Float64}}
    edge_list::Vector{Tuple{Int64,Int64}}
    cell_area::Vector{Float64}
    cell_vertices::Vector{Vector{Tuple{Float64,Float64}}}
    border::Vector{Tuple{Float64,Float64}}
    scale_factor::Float64
    Nedge::Int64
    Nnode::Int64
end

mutable struct ContModel
    mesh::Mesh
    id_slack::Int64
    p::Vector{Float64}
    dp
    b::Vector{Float64}
    m::Vector{Float64}
    d::Vector{Float64}
    th::Vector{Float64}
end

include("disturbances.jl")
include("dynamics.jl")
include("mesh.jl")
include("param.jl")
include("plot.jl")
include("stable.jl")
include("utils.jl")

end
