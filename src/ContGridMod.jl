module ContGridMod

using BlockArrays
using DifferentialEquations
using Ferrite
using FerriteGmsh
using FerriteViz
using FileIO
using CairoMakie
using Gmsh
using HDF5
using IterativeSolvers
using JLD2
using JSON
using LinearAlgebra
using Plots
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

mutable struct ContModel
    grid::Grid
    dh₁::DofHandler
    dh₂::DofHandler
    cellvalues::CellScalarValues
    area::Float64
    m_nodal::Vector{Float64}
    d_nodal::Vector{Float64}
    p_nodal::Vector{Float64}
    bx_nodal::Vector{Float64}
    by_nodal::Vector{Float64}
    θ₀_nodal::Vector{Float64}
    fault_nodal::Vector{Float64}
    m::Function
    d::Function
    p::Function
    bx::Function
    by::Function
    θ₀::Function
    fault::Function
    ch::ConstraintHandler
    K₀::SparseMatrixCSC
    f₀::Vector{Float64}
end

include("discrete.jl");
include("disturbances.jl")
include("dynamics.jl")
include("mesh.jl")
include("params.jl")
include("plotting.jl")
include("stable.jl")
include("tools.jl")

end
