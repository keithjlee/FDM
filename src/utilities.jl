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

    DataInspector(fig)
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