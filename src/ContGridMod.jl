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
using JSON
using LinearAlgebra
using Plots
using SparseArrays


abstract type GridModel end

mutable struct DiscModel <: GridModel
    m_gen::Vector{<:Real}
    d_gen::Vector{<:Real}
    id_gen::Vector{<:Int}
    id_slack::Int
    coord::Matrix{<:Real}
    d_load::Vector{<:Real}
    id_line::Matrix{Int}
    b::Vector{<:Real}
    p_load::Vector{<:Real}
    th::Vector{<:Real}
    p_gen::Vector{<:Real}
    max_gen::Vector{<:Real}
    Nbus::Int
    Ngen::Int
    Nline::Int
end

mutable struct ContModel <: GridModel
    grid::Grid
    dh₁::DofHandler
    dh₂::DofHandler
    cellvalues::CellScalarValues
    area::Real
    m_nodal::Vector{<:Real}
    d_nodal::Vector{<:Real}
    p_nodal::Vector{<:Real}
    bx_nodal::Vector{<:Real}
    by_nodal::Vector{<:Real}
    θ₀_nodal::Vector{<:Real}
    fault_nodal::Vector{<:Real}
    m::Function
    d::Function
    p::Function
    bx::Function
    by::Function
    θ₀::Function
    fault::Function
    ch::ConstraintHandler
    ContModel(
        grid::Grid,
        dh₁::DofHandler,
        dh₂::DofHandler,
        cellvalues::CellScalarValues,
        area::Real,
        m_nodal::Vector{<:Real},
        d_nodal::Vector{<:Real},
        p_nodal::Vector{<:Real},
        bx_nodal::Vector{<:Real},
        by_nodal::Vector{<:Real},
        θ₀_nodal::Vector{<:Real},
        fault_nodal::Vector{<:Real},
        ch::ConstraintHandler) = new(
        grid,
        dh₁,
        dh₂,
        cellvalues,
        area,
        m_nodal,
        d_nodal,
        p_nodal,
        bx_nodal,
        by_nodal,
        θ₀_nodal,
        fault_nodal,
        (x; extrapolate=true, warn=:semi) -> interpolate(x, grid, dh₁, m_nodal, :u, extrapolate=extrapolate, warn=warn),
        (x; extrapolate=true, warn=:semi) -> interpolate(x, grid, dh₁, d_nodal, :u, extrapolate=extrapolate, warn=warn),
        (x; extrapolate=true, warn=:semi) -> interpolate(x, grid, dh₁, p_nodal, :u, extrapolate=extrapolate, warn=warn),
        (x; extrapolate=true, warn=:semi) -> interpolate(x, grid, dh₁, bx_nodal, :u, extrapolate=extrapolate, warn=warn),
        (x; extrapolate=true, warn=:semi) -> interpolate(x, grid, dh₁, by_nodal, :u, extrapolate=extrapolate, warn=warn),
        (x; extrapolate=true, warn=:semi) -> interpolate(x, grid, dh₁, θ₀_nodal, :u, extrapolate=extrapolate, warn=warn),
        (x; extrapolate=true, warn=:semi) -> interpolate(x, grid, dh₁, fault_nodal, :u, extrapolate=extrapolate, warn=warn),
        ch
    )
    ContModel(; grid::Grid,
        dh₁::DofHandler,
        dh₂::DofHandler,
        cellvalues::CellScalarValues,
        area::Real,
        m_nodal::Vector{<:Real},
        d_nodal::Vector{<:Real},
        p_nodal::Vector{<:Real},
        bx_nodal::Vector{<:Real},
        by_nodal::Vector{<:Real},
        θ₀_nodal::Vector{<:Real},
        fault_nodal::Vector{<:Real},
        ch::ConstraintHandler) = ContModel(grid,
        dh₁::DofHandler,
        dh₂,
        cellvalues,
        area,
        m_nodal,
        d_nodal,
        p_nodal,
        bx_nodal,
        by_nodal,
        θ₀_nodal,
        fault_nodal,
        ch)
end

include("discrete.jl")
include("disturbances.jl")
include("dynamics.jl")
include("mesh.jl")
include("params.jl")
include("plotting.jl")
include("stable.jl")
include("tools.jl")

end
