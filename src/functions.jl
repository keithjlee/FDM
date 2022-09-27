using LinearAlgebra, kjlMakie, GLMakie, Statistics

set_theme!(kjl_light)

"""
Extracts the vector of Boolean DOFs in node sequence.
|dofs| = |nodes|
"""
function dofs(points::Vector{Node})
    return [p.dof for p in points]
end

"""
Extracts the indices of free (N) and fixed (F) nodes
"""
function NF(d::Vector{Bool})
    return findall(d), findall(.!d)
end

"""
Extracts the fixed/free indices of a vector of nodes
"""
function NF(points::Vector{Node})
    d = dofs(points)
    return NF(d)
end

"""
Directly on network
"""
function NF!(network::Network)
    network.N, network.F = NF(network.nodes)
end

"""
Creates a 1x3 vector of nodal position
"""
function deconstructNode(point::Node; horizontal = true)
    if horizontal
        return Float64.([point.x point.y point.z])
    else
        return Float64.([point.x, point.y, point.z])
    end
end

"""
Creates an nx3 matrix of nodal positions
"""
function deconstructNode(points::Vector{Node})
    x = [p.x for p in points]
    y = [p.y for p in points]
    z = [p.z for p in points]

    return Float64.([x y z])
end

"""
Directly on a network
"""
function deconstructNodes!(network::Network)
    network.xyz = deconstructNode(network.nodes)
end


"""
Creates load matrix P (nx3) with loads indexed vertically w/r/t node order
"""
function loadMatrix(loads::Vector{Load}, n::Int64)
    p = zeros(n, 3)
    for load in loads
        p[load.index, :] = load.force
    end
    return p
end

"""
directly on a network
"""
function loadMatrix!(network::Network)
    network.P = loadMatrix(network.loads, length(network.nodes))
end



"""
Preprocessing of data for a new FDM network
"""
function process!(network::Network)
    #fixed-free indices N, F
    NF!(network)

    #branch node matrix C
    branchMatrix!(network)
    network.Cn = network.C[:, network.N]
    network.Cf = network.C[:, network.F]

    #nodal positions [x y z]
    deconstructNodes!(network)

    #force density vector q
    forceDensities!(network)
    network.Q = diagm(network.q)

    #load matrix P
    loadMatrix!(network)
    network.Pn = network.P[network.N, :]

    network.processed = true
end

"""
performs the fdm analysis
"""
function solve!(network::Network; reprocess = false)
    if !network.processed || reprocess
        process!(network)
    end

    # solve for final free node positions
    network.xyz[network.N, :] = (network.Cn' * network.Q * network.Cn) \ (network.Pn - network.Cn' * network.Q * network.Cf * network.xyz[network.F, :])
end

"""
Extracts force density vector (q) from elements
"""
function forceDensities(elements::Vector{Element})
    return [e.q for e in elements]
end

function forceDensities!(network::Network)
    network.q = forceDensities(network.elements)
end

"""
Creates the branch-node connectivity matrix C
"""
function branchMatrix(elements::Vector{Element}, points::Vector{Node})

    #initialize
    c = zeros(Int64, length(elements), length(points))

    #rows = elements, columns = nodes
    for (i, element) in enumerate(elements)
        c[i, element.iStart] = -1
        c[i, element.iEnd] = 1
    end

    return c
end

"""
Directly on a network
"""
function branchMatrix!(network::Network)
    network.C = branchMatrix(network.elements, network.nodes)
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
Plot a network (must be at least processed, but solved is optional)
"""
function plot(network::Network)
    if !network.processed
        error("run process!(network) first.")
    end

    # initialize
    fig = Figure(backgroundcolor = :white)
    ax = Axis3(fig[1,1],
        aspect = :data)

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
