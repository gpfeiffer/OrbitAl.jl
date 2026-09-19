module OrbitAl

include("permutation.jl")
using .permutation
export Perm, degree, domain, cycles, shape, order, isidentity, last_moved
export permuted, transpositions, shuffle!

include("orbits.jl")
using .orbits
export Orbit
export orbit, orbitl, onPoints, onRight, onWords, onPairs, onSets
export orbit_with_words, orbit_with_transversal, orbit_with_stabilizer
export orbit_with_dist, orbit_with_tree, orbit_with_edges, orbit_with_images
export orbitx, orbitx_with_words, orbitx_with_dist, orbitx_with_edges
export edges_from_images
export orbit_and_more, orbit_and_data
export Item, orbit_with_data

include("permgroup.jl")
using .permgroup
export PermGp, elements, conjClasses, closure, subgroups, subgpClasses
export sizeOfGroup, randomGroupElement, memberOfGroup
export isPrimePower, zuppos

include("bfsdfs.jl")
include("plotting.jl")

##  Methods are added by the extensions in ext/, when the packages they need
##  are loaded.
plot_edges(args...) =
    error("plot_edges needs Graphs and GraphPlot: `using Graphs, GraphPlot`")
write_d3_edges(args...) = error("write_d3_edges needs JSON: `using JSON`")
write_d3_col_edges(args...) = error("write_d3_col_edges needs JSON: `using JSON`")

export plot_edges, write_d3_edges, write_d3_col_edges

# on-demand modules: load with e.g. `using OrbitAl.syt`
include("coxeter.jl")
include("syt.jl")
include("simsgroup.jl")
include("shifts.jl")
include("involution.jl")
include("coset.jl")
include("presentations.jl")
include("variants.jl")
include("enumerator.jl")

end # module OrbitAl
