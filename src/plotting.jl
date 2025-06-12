#############################################################################
##
#A  plotting.jl                                                       OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  Plotting graphs in jupyter notebooks using different interfaces.
##
module plotting

using Graphs, GraphPlot

export plot_edges

function plot_edges(edges)
    graph = SimpleGraph(Edge.(edges))
    gplot(graph, nodelabel=vertices(graph))
end


end # module
