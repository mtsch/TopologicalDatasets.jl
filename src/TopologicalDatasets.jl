module TopologicalDatasets

using Base.Cartesian
using Distances
using LightGraphs
using LinearAlgebra
using MultivariateStats
using NearestNeighbors
using Random
using RecipesBase
using SimpleWeightedGraphs
using SparseArrays
using StaticArrays
using StatsBase

using Base: @kwdef

export GeodesicMetric, distances, points
export Ball, Cube, Klein, Knot, Noisy, Sphere, Torus, AsymTorus, sample, ×

include("generators.jl")
include("geodesicmetrics.jl")

end
