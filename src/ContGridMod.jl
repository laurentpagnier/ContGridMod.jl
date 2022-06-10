module ContGridMod

using SparseArrays

mutable struct DiscModel
    mg::Array{Float64, 1}
    dg::Array{Float64, 1}
    idgen::Array{Int64, 1}
    coord::Array{Float64, 2}
    dl::Array{Float64, 1}
    idb::Array{Int64, 2}
    bline::Array{Float64, 1}
    load::Array{Float64, 1}
    th::Array{Float64, 1}
    gen::Array{Float64, 1}
    max_gen::Array{Float64, 1}
    Nbus::Int64
end

mutable struct Mesh
    Nx::Int64
    Ny::Int64
    coord::Matrix{Float64}
    line_coord::Matrix{Float64}
    inc_mat::Matrix{Int64}
    isgrid::BitVector
    yrange::Vector{Float64}
    xrange::Vector{Float64}
    dx::Float64
end

mutable struct ContModel
    mesh::Mesh
    minv::Vector{Float64}
    gamma::Vector{Float64}
    p::Vector{Float64}
    xi::SparseMatrixCSC{Float64, Int64}
    b::Vector{Float64}
    m::Vector{Float64}
    d::Vector{Float64}
    th::Vector{Float64}
    scale_factor::Float64
    dmax::Float64
    Niter::Int64
    tau::Float64
    patch::Float64
    min_factor::Float64
end

include("disc_solvers.jl");
include("disturbances.jl")
include("dynamics.jl")
include("mesh.jl")
include("params.jl")
include("plotting.jl")
include("ps_analysis.jl")
include("stable.jl")
include("tools.jl")

end
