"""
Node in fdm network. Includes positional data and fixed/free information
"""
mutable struct Node
    x::Union{Int64, Float64}
    y::Union{Int64, Float64}
    z::Union{Int64, Float64}
    dof::Bool # true = free; false = fixed
    id::Union{Symbol, Nothing}

    #empty constructor
    function Node()
        return new()
    end

    # individual coordinate basis
    function Node(x::Union{Int64, Float64}, y::Union{Int64, Float64}, z::Union{Int64, Float64}, dof::Bool)
        return new(x, y, z, dof, nothing)
    end

    # using a vector to represent position
    function Node(pos::Union{Vector{Int64}, Vector{Float64}}, dof::Bool)
        if length(pos) != 3
            error("pos should be length 3.")
        else
            return new(pos..., dof, nothing)
        end
    end

    #with an id
    function Node(x::Union{Int64, Float64}, y::Union{Int64, Float64}, z::Union{Int64, Float64}, dof::Bool, id::Symbol)
        return new(x, y, z, dof, id)
    end

    function Node(pos::Union{Vector{Int64}, Vector{Float64}}, dof::Bool, id::Symbol)
        if length(pos) != 3
            error("pos should be length 3.")
        else
            return new(pos..., dof, id)
        end
    end
end

"""
Element represents an edge in the network. Includes the global indices and positional information of the start and end nodes, as well as the force density.
"""
mutable struct Element
    pStart::Node #start point
    iStart::Int64 #index of start point in vector of points
    pEnd::Node #end point
    iEnd::Int64 #index of end point in vector of points
    q::Union{Int64, Float64} #force density
    id::Union{Symbol, Nothing}

    function Element(points::Vector{Node}, iStart::Int64, iEnd::Int64, q::Union{Int64, Float64})
        element = new(points[iStart], iStart, points[iEnd], iEnd, q, nothing)
        return element
    end

    function Element(points::Vector{Node}, iStart::Int64, iEnd::Int64, q::Union{Int64, Float64}, id::Symbol)
        element = new(points[iStart], iStart, points[iEnd], iEnd, q, id)
        return element
    end
end

"""
Load type is associated with a node (and its position in node vector), as well as a length 3 force vector (x,y,z)
"""
mutable struct Load
    point::Node # point at which load is applied
    index::Int64 # position of point in vector of points
    force::Union{Vector{Int64}, Vector{Float64}} # force vector

    function Load(points::Vector{Node}, i::Int64, force::Union{Vector{Int64}, Vector{Float64}})
        if length(force) != 3
            error("Force vector should be length 3.")
        else
            load = new(points[i], i, force)
            return load
        end
    end
end

"""
An FDM network with all relevant information for analysis
"""
mutable struct Network
    nodes::Vector{Node} #vector of nodes
    elements::Vector{Element} #vector of elements
    loads::Vector{Load} #vector of loads
    q::Union{Vector{Float64}, Vector{Int64}} #vector of force densities
    Q::Union{SparseMatrixCSC{Float64, Int64}, SparseMatrixCSC{Int64, Int64}} #diagm(q)
    C::SparseMatrixCSC{Int64, Int64} #branch node matrix
    N::Vector{Int64} #fixed node indices
    F::Vector{Int64} #free node indices
    Cn::SparseMatrixCSC{Int64, Int64} #fixed branch node matrix
    Cf::SparseMatrixCSC{Int64, Int64} #free branch node matrix
    P::Union{Matrix{Float64}, Matrix{Int64}} #load matrix
    Pn::Union{Matrix{Float64}, Matrix{Int64}} #free node load matrix
    xyz::Union{Matrix{Int64}, Matrix{Float64}} #xyz matrix of initial positions
    # xyzn::Union{Matrix{Int64}, Matrix{Float64}} #xyz matrix of new positions
    # xyzf::Union{Matrix{Int64}, Matrix{Float64}} # xyz matrix of fixed positions
    processed::Bool

    function Network(nodes::Vector{Node}, elements::Vector{Element}, loads::Vector{Load})
        network = new(nodes, elements, loads)
        network.processed = false
        return network
    end
end
