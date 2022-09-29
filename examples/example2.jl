using FDM

# create grid of nodes
begin
    xrange = collect(0:20);
    yrange = collect(0:20);

    nodes = [Node(x, y, 0, true) for x in xrange, y in yrange];

    for node in nodes
        if node.x == 0 && (node.y == 0 || node.y == 20)
            node.dof = false
        end

        if node.x == 20 && (node.y == 0 || node.y == 20)
            node.dof = false
        end
    end
end

# create the grid of elements (this is a bit hacky)
begin
    nodeMatIndex = LinearIndices(nodes);
    i, j = size(nodes);
    nodes = vcat(nodes...);

    elements = Vector{Element}();

    for x = 1:i-1
        for y = 1:j-1
            push!(elements, Element(nodes, nodeMatIndex[x,y], nodeMatIndex[x,y+1], 10.))
            push!(elements, Element(nodes, nodeMatIndex[x,y], nodeMatIndex[x+1,y], 10.))
        end
    end;

    for x = 1:j-1
        push!(elements, Element(nodes, nodeMatIndex[x, j], nodeMatIndex[x+1, j], 10.))
        push!(elements, Element(nodes, nodeMatIndex[j, x], nodeMatIndex[j, x+1], 10.))
    end;
end

#create loads (only applied to free nodes)
begin
    loads = Vector{Load}();

    for (x, node) in enumerate(nodes)
        if node.dof == true
            push!(loads, Load(nodes, x, [0., 0., -1.]))
        end
    end;
end

#assemble network
network = Network(nodes, elements, loads);

#preprocess network
process!(network);

#analyze for nodal positions
solve!(network);

#update nodal positions to network.nodes
xyzUpdate!(network);

#visualize
networkFig = plot(network)