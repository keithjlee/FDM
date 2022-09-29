module FDM

using LinearAlgebra
using kjlMakie
using GLMakie
using Statistics

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
export plot

end # module FDM
