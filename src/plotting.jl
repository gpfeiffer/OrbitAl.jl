#############################################################################
##
#A  plotting.jl                                                       OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  Plotting graphs in jupyter notebooks using different interfaces.
##
module plotting

export plot_edges, write_d3_edges, write_d3_col_edges

using Graphs, GraphPlot

function plot_edges(edges)
    graph = SimpleGraph(Edge.(edges))
    gplot(graph, nodelabel=vertices(graph))
end

# using Compose, Cairo
# p = plot_edges(orb.edges)
# draw(PDF("graph.pdf",  600px, 400px), p)

using JSON

function d3_json(edges)
    nodes = Set(vcat(collect.(edges)...))
    links = [Dict("source" => string(i), "target" => string(j)) for (i, j) in edges]
    Dict("nodes" => [Dict("id" => string(n)) for n in nodes], "links" => links)
end

function d3_col_json(edges)
    nodes = union([Set(e[[1,2]]) for e in edges]...)
    links = [
        Dict("source" => i, "target" => j, "label" => k)
        for (i, j, k) in edges
    ]
    Dict("nodes" => [Dict("id" => n) for n in nodes], "links" => links)
end

using Base64, IJulia

function html_d3_force_graph(graph_json::String)
    return """
    <div id="d3graph"></div>
    <script src="https://d3js.org/d3.v7.min.js"></script>
    <script>
    const width = 600, height = 400;

    const graph = $graph_json;

    const svg = d3.select("#d3graph")
        .append("svg")
        .attr("width", width)
        .attr("height", height);

    const simulation = d3.forceSimulation(graph.nodes)
        .force("link", d3.forceLink(graph.links).id(d => d.id).distance(50))
        .force("charge", d3.forceManyBody().strength(-200))
        .force("center", d3.forceCenter(width / 2, height / 2));

    const link = svg.append("g")
        .attr("stroke", "#aaa")
        .selectAll("line")
        .data(graph.links)
        .join("line");

    const node = svg.append("g")
        .attr("stroke", "#fff")
        .attr("stroke-width", 1.5)
        .selectAll("circle")
        .data(graph.nodes)
        .join("circle")
        .attr("r", 8)
        .attr("fill", "steelblue")
        .call(drag(simulation));

    const label = svg.append("g")
        .selectAll("text")
        .data(graph.nodes)
        .join("text")
        .text(d => d.id)
        .attr("font-size", 12)
        .attr("dy", -10);

    simulation.on("tick", () => {
        link
            .attr("x1", d => d.source.x)
            .attr("y1", d => d.source.y)
            .attr("x2", d => d.target.x)
            .attr("y2", d => d.target.y);

        node
            .attr("cx", d => d.x)
            .attr("cy", d => d.y);

        label
            .attr("x", d => d.x)
            .attr("y", d => d.y);
    });

    function drag(simulation) {
        function dragstarted(event, d) {
            if (!event.active) simulation.alphaTarget(0.3).restart();
            d.fx = d.x;
            d.fy = d.y;
        }

        function dragged(event, d) {
            d.fx = event.x;
            d.fy = event.y;
        }

        function dragended(event, d) {
            if (!event.active) simulation.alphaTarget(0);
            d.fx = null;
            d.fy = null;
        }

        return d3.drag()
            .on("start", dragstarted)
            .on("drag", dragged)
            .on("end", dragended);
    }
    </script>
    """
end

function html_d3_col_force_graph(graph_json::String)
    return """
    <div id="d3graph"></div>
    <script src="https://d3js.org/d3.v7.min.js"></script>
    <script>
    const width = 1200, height = 800;

    const graph = $graph_json;

    const color = d3.scaleOrdinal(d3.schemeCategory10);

    const svg = d3.select("#d3graph")
        .append("svg")
        .attr("width", width)
        .attr("height", height);

    const simulation = d3.forceSimulation(graph.nodes)
        .force("link", d3.forceLink(graph.links).id(d => d.id).distance(5))
        .force("charge", d3.forceManyBody().strength(-40))
        .force("center", d3.forceCenter(width / 2, height / 2));

    const link = svg.append("g")
        .attr("stroke", "#aaa")
        .selectAll("line")
        .data(graph.links)
        .enter().append("line")
        .attr("class", "link")
        .attr("stroke", d => color(d.label))
        .attr("stroke-width", 1);

    const node = svg.append("g")
        .attr("stroke", "#fff")
        .attr("stroke-width", 1.5)
        .selectAll("circle")
        .data(graph.nodes)
        .join("circle")
        .attr("r", 5)
        .attr("fill", "steelblue")
        .call(drag(simulation));

    const label = svg.append("g")
        .selectAll("text")
        .data(graph.nodes)
        .join("text")
        .text(d => d.id)
        .attr("font-size", 12)
        .attr("dy", -10);

    simulation.on("tick", () => {
        link
            .attr("x1", d => d.source.x)
            .attr("y1", d => d.source.y)
            .attr("x2", d => d.target.x)
            .attr("y2", d => d.target.y);

        node
            .attr("cx", d => d.x)
            .attr("cy", d => d.y);

        label
            .attr("x", d => d.x)
            .attr("y", d => d.y);
    });

    function drag(simulation) {
        function dragstarted(event, d) {
            if (!event.active) simulation.alphaTarget(0.3).restart();
            d.fx = d.x;
            d.fy = d.y;
        }

        function dragged(event, d) {
            d.fx = event.x;
            d.fy = event.y;
        }

        function dragended(event, d) {
            if (!event.active) simulation.alphaTarget(0);
            d.fx = null;
            d.fy = null;
        }

        return d3.drag()
            .on("start", dragstarted)
            .on("drag", dragged)
            .on("end", dragended);
    }
    </script>
    """
end

function write_d3_edges(edges, filename::String="graph.html")
    graph_data = d3_json(edges)
    graph_json = JSON.json(graph_data)
    html = html_d3_force_graph(graph_json)
    write(filename, html)
end

function write_d3_col_edges(edges, filename::String="graph.html")
    graph_data = d3_col_json(edges)
    graph_json = JSON.json(graph_data)
    html = html_d3_col_force_graph(graph_json)
    write(filename, html)
end

end # module
