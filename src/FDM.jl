module FDM

using LinearAlgebra
using kjlMakie
using GLMakie
using Statistics
using SparseArrays
using LinearSolve

include("types.jl")
export Node
export Element
export Load
export Network

include("functions.jl")
include("utilities.jl")
export deconstructNode
export process!
export solve!
export forceDensities
export dist
export memberForces
export memberLengths
export FL
export plot
export qUpdate!
export xyzUpdate!

end # module FDM
