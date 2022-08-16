module ContGridMod

using JSON
using SparseArrays

mutable struct DiscModel
    m_gen::Vector{Float64} # generator inertia
    d_gen::Vector{Float64} # generator damping
    id_gen::Vector{Int64} # generator inertia
    id_slack::Int64 # index of the slack bus
    coord::Matrix{Float64} # coordinates in 
    d_load::Vector{Float64}
    id_line::Matrix{Int64} # list of indices
    b::Vector{Float64} # susceptance
    p_load::Vector{Float64} 
    th::Vector{Float64}
    p_gen::Vector{Float64}
    max_gen::Vector{Float64}
    Nbus::Int64
    Ngen::Int64
    Nline::Int64
end

mutable struct Mesh
    Nx::Int64
    Ny::Int64
    coord::Matrix{Float64}
    line_coord::Matrix{Float64}
    id_edge::Matrix{Int64}
    is_grid::BitVector
    yrange::Vector{Float64}
    xrange::Vector{Float64}
    dx::Float64
    border::Matrix{Float64}
    scale_factor::Float64
    Nedge::Int64
    Nnode::Int64
end

mutable struct ContModel
    mesh::Mesh
    id_slack::Int64
    p::Vector{Float64}
    dp::Vector{Float64}
    b::Vector{Float64}
    m::Vector{Float64}
    d::Vector{Float64}
    th::Vector{Float64}
end

include("discrete.jl");
include("disturbances.jl")
include("dynamics.jl")
include("mesh.jl")
include("params.jl")
include("plotting.jl")
include("ps_analysis.jl")
include("stable.jl")
include("tools.jl")

end
