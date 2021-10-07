module ContGridMod

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

mutable struct ContModel
    Nx::Int64
    Ny::Int64
    coord::Array{Float64, 2}
    isinside::BitVector
    isborder::BitVector
    isgrid::BitVector
    yrange::Array{Float64, 1}
    xrange::Array{Float64, 1}
    n::Array{Float64, 2}
    dx::Float64
    minv::Array{Float64, 1}
    gamma::Array{Float64, 1}
    p::Array{Float64, 1}
    xi::SparseMatrixCSC{Float64, Int64}
    bx::Array{Float64, 1}
    by::Array{Float64, 1}
    m::Array{Float64, 1}
    d::Array{Float64, 1}
    th::Array{Float64, 1}
end

mutable struct Mesh
    Nx::Int64
    Ny::Int64
    coord::Array{Float64, 2}
    isinside::BitVector
    isborder::BitVector
    isgrid::BitVector
    yrange::Array{Float64, 1}
    xrange::Array{Float64, 1}
    n::Array{Float64, 2}
    dx::Float64
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
