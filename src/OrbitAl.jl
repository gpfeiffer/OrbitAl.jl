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
export intersect_groups

include("bfsdfs.jl")
include("plotting.jl")
include("syt.jl")
using .syt
export newtonDif, newtonSum, newtonDifR, newtonSumR
export compositionSubset, subsetComposition
export takeAway, subsets, partitions, standardYTs, tableau_path
include("simsgroup.jl")

include("coxeter.jl")
using .coxeter
export coxeterGraph, coxeterMat, cartanMat, CoxeterGp, coxeterConjugacyClasses
export coxeterLength, coxeterWord, permCoxeterWord, reflections
export prefixes, prefixes_with_edges, shapes, longestElt, parabolicTransversal

include("shifts.jl")
include("involution.jl")
using .involution
export onInvolutions, involutions, onInvolutionClasses, involutionClasses

include("presentations.jl")
include("variants.jl")
include("enumerator.jl")

end # module OrbitAl
