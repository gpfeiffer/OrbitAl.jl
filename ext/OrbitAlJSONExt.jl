#############################################################################
##
#A  OrbitAlJSONExt.jl                                                 OrbitAl
##
#C  The d3 writers, when JSON is loaded.
##
module OrbitAlJSONExt

using OrbitAl, JSON

function OrbitAl.write_d3_edges(edges, filename::String="graph.html")
    json = JSON.json(OrbitAl.plotting.d3_json(edges))
    write(filename, OrbitAl.plotting.html_d3_force_graph(json))
end

function OrbitAl.write_d3_col_edges(edges, filename::String="graph.html")
    json = JSON.json(OrbitAl.plotting.d3_col_json(edges))
    write(filename, OrbitAl.plotting.html_d3_col_force_graph(json))
end

end # module OrbitAlJSONExt
