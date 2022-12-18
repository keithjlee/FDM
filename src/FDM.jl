module FDM

using LinearAlgebra
using kjlMakie
using GLMakie
using Statistics
using SparseArrays

include("types.jl")
export Node
export Element
export Load
export Network

using JSON
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
export initialLengths
export FDMremoteread
export FDMremotesave

end # module FDM
