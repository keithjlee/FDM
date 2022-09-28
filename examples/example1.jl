using FDM

# define Nodes
begin
    P0 = Node([0, 0, 0], false)
    P1 = Node([1, -2, 0], true)
    P2 = Node([2, -1, 0], true)
    P3 = Node([3, 0.5, 0], false)
    P4 = Node([5.2, -1, 0], false)
    P5 = Node([3, -2, 0], true)
    P6 = Node([3, -3, 0], false)

    nodes = [P0,
        P1,
        P2,
        P3,
        P4,
        P5,
        P6]
end

# define elements
begin
    indices = [(0,3),
        (3, 4),
        (2,3),
        (0,2),
        (0,1),
        (2,5),
        (1,6),
        (4,5),
        (5,6),
        (1,2)]

    elements = Vector{Element}()
    for (i, index) in enumerate(indices)
        push!(elements, Element(nodes, index[1] + 1, index[2] + 1, i))
    end

end

#define loads
loads = Vector{Load}()
for (i, node) in enumerate(nodes)
    if node.dof
        push!(loads, Load(nodes, i , [0., 0., 10.]))
    end
end

#create network
network = Network(nodes, elements, loads)

#preprocess (optional)
process!(network)

#solve
solve!(network)

#visualize
networkFig = plot(network)