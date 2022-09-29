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
    network.Q = spdiagm(network.q)

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
    c = spzeros(Int64, length(elements), length(points))

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

