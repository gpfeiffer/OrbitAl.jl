module OrbitAl

include("permutation.jl")
using .permutation
export Perm, degree, domain, cycles, shape, order, isidentity, last_moved
export permuted, transposition, transpositions, shuffle!

include("orbits.jl")
using .orbits

include("bfsdfs.jl")
include("syt.jl")

end # module OrbitAl
