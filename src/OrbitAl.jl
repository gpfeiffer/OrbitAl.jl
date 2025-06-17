module OrbitAl

include("permutation.jl")
using .permutation
export Perm, degree, domain, cycles, shape, order, isidentity, last_moved
export permuted, transposition, transpositions, shuffle!

include("orbits.jl")
using .orbits

include("permgroup.jl")
using .permgroup
export PermGp, elements, conjClasses, closure, subgroups, subgpClasses
export sizeOfGroup, randomGroupElement, memberOfGroup

include("bfsdfs.jl")
include("plotting.jl")
include("syt.jl")

end # module OrbitAl
