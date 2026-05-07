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

include("permgroup.jl")
using .permgroup
export PermGp, elements, conjClasses, closure, subgroups, subgpClasses
export sizeOfGroup, randomGroupElement, memberOfGroup

include("bfsdfs.jl")
include("plotting.jl")
include("syt.jl")
include("simsgroup.jl")

include("coxeter.jl")
using .coxeter
export coxeterGraph, coxeterMat, cartanMat, CoxeterGp, coxeterConjugacyClasses
export coxeterLength, coxeterWord, permCoxeterWord, reflections
export prefixes, prefixes_with_edges, shapes, longestElt, parabolicTransversal

include("shifts.jl")
include("involution.jl")

include("presentations.jl")
include("variants.jl")
include("enumerator.jl")

end # module OrbitAl
