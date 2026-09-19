#############################################################################
##
#A  OrbitAlGraphPlotExt.jl                                            OrbitAl
##
#C  plot_edges, when Graphs and GraphPlot are loaded.
##
module OrbitAlGraphPlotExt

using OrbitAl, Graphs, GraphPlot

function OrbitAl.plot_edges(edges)
    graph = SimpleGraph(Edge.(edges))
    gplot(graph, nodelabel=vertices(graph))
end

# using Compose, Cairo
# p = plot_edges(orb.edges)
# draw(PDF("graph.pdf",  600px, 400px), p)

end # module OrbitAlGraphPlotExt
