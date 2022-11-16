"""
Plot a network (must be at least processed, but solved is optional)
"""
function plot(network::Network; showaxis = true)
    if !network.processed
        error("run process!(network) first.")
    end

    # initialize
    fig = Figure(backgroundcolor = :white)
    ax = Axis3(fig[1,1],
        aspect = :data)
    
    if !showaxis
        hidedecorations!(ax)
        hidespines!(ax)
    end

    # interactivity
    sg = SliderGrid(fig[2,1],
        (label = "Line Scale", range = 0.1:0.1:10, startvalue = 1.),
        (label = "Point Scale", range = 0.1:0.1:10, startvalue = 1.))

    q = forceDensities(network.elements)

    # normalize
    q ./= maximum(q)

    linefactor = lift(sg.sliders[1].value) do v
        v .* q
    end

    nodefactor = lift(sg.sliders[2].value) do v
        0.1 * v
    end

    # Nodes as primitives
    points = Point3.(eachrow(network.xyz))

    # elements
    l = vcat([[points[e.iStart], points[e.iEnd]] for e in network.elements]...)
    lines = linesegments!(l, 
        color = (:black, 0.75),
        linewidth = linefactor)

    # node colors: black = fixed; blue = free
    colors = [node.dof ? blue : :black for node in network.nodes]

    joints = meshscatter!(points,
        color = colors,
        markersize = nodefactor)

    # loads
    arrs = Vec3.(eachrow(network.P))
    loadLengths = norm.(arrs)
    elemLengths = norm.(eachrow(network.C * network.xyz))
    scaleFactor = mean(elemLengths) / maximum(loadLengths)

    forces = arrows!(points, arrs .* scaleFactor,
        color = pink)

    # options
    toggles = [Toggle(fig, active = true)]
    labels = [Label(fig, "Loads")]

    fig[1,2] = grid!(hcat(toggles, labels), tellheight = false)
    connect!(forces.visible, toggles[1].active)

    # DataInspector(fig)
    # show figure
    display(fig)


    return fig
end

"""
Euclidean distance metric for Node types
"""
function dist(a::Node, b::Node)
    return sqrt((a.x - b.x)^2 + (a.y - b.y)^2 + (a.z - b.z)^2)
end

"""
Extract member forces
"""
function memberForces(network::Network)
    return norm.(eachrow(network.C * network.xyz)) .* network.q
end

"""
Extract member lengths
"""
function memberLengths(network::Network)
    return norm.(eachrow(network.C * network.xyz))
end

"""
Custom indexing based on IDs of structs
"""
function Base.getindex(nodes::Vector{Node}, i::Symbol)
    return [node for node in nodes if node.id == i]
end

"""
Custom indexing for elements
"""
function Base.getindex(elements::Vector{Element}, i::Symbol)
    return [element for element in elements if element.id == i]
end

"""
Custom indexing of networks
"""
function Base.getindex(network::Network, i::Symbol)
    nodes = network.nodes[i]
    elements = network.elements[i]

    return nodes, elements
end

"""
findall methods for querying
"""
function Base.findall(elements::Vector{Element}, i::Symbol)
    return findall([x.id == i for x in elements])
end

function Base.findall(nodes::Vector{Node}, i::Symbol)
    return findall([x.id == i for x in nodes])
end

function Base.findall(network::Network, i::Symbol)
    return findall(network.elements, i), findall(network.nodes, i)
end

"""
Repopulate network with new q values
"""
function qUpdate!(network::Network, q::Union{Vector{Float64}, Vector{Int64}})
    if length(network.elements) != length(q)
        error("Number of elements and q must be equal.")
    end

    for (i, element) in enumerate(network.elements)
        element.q = q[i]
    end

    forceDensities!(network)
end

#single value
function qUpdate!(network::Network, q::Union{Float64, Int64})
    for element in network.elements
        element.q = q
    end

    forceDensities!(network)
end

#by id
function qUpdate!(network::Network, q::Union{Float64, Int64}, id::Symbol)
    for element in network.elements[id]
        element.q = q
    end

    forceDensities!(network)
end

"""
Get initial lengths for form finding;
assumes units of E, A are consistent with L
"""
function initialLengths(network::Network, E::Union{Float64, Int64}, A::Union{Float64, Int64})
    n = length(network.elements) #number of elements
    Id = I(n) #identity matrix n × n
    Em = E * Id #diagonal matrix of stiffness, E
    Am = A * Id #diagonal matrix of areas, A
    L = spdiagm(memberLengths(network)) #diagonal matrix of final lengths

    return diag((Id + (Em * Am) \ network.Q * L) \ Id)
end

"""
Initial length method for varying section properties
"""
function initialLengths(network::Network, E::Union{Vector{Float64}, Vector{Int64}}, A::Union{Vector{Float64}, Vector{Int64}})
    n = length(network.elements)

    # make sure material property vectors are the same
    @assert n == length(E) == length(A) "E and A vectors must be equal length"

    # matrix representation of components
    Id = I(n)
    Em = diagm(E)
    Am = diagm(A)
    L = spdiagm(memberLengths(network))

    return diag((Id + (Em * Am) \ network.Q * L) \ Id)
end