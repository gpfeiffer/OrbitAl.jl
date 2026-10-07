#############################################################################
##
#A  linear.jl                                                         OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  Linear orbit algorithms: spinning
##
module linear

using ..sparsevec
using ..unionfind

export spinning

"""
    spinning(aaa, x, under)

The spinning algorithm, the linear version of the orbit algorithm: returns a
basis of the smallest subspace that contains the vector `x` and is invariant
under the operators in `aaa`, acting via `under`.  The basis consists of `x`
and images of earlier basis vectors.

A vector is new if it is not in the span of the basis so far.  Linear
Union-Find tests this incrementally: each basis vector `v` is entered as the
relation `v = 0`, and `unite!` returns `0` exactly for the vectors in the span.
"""
function spinning(aaa, x, under)

    # the span of list, as relations v = 0
    parent = Dict{Int, SparseVec{eltype(x)}}()

    # records v, if new
    function inSpan!(v)
        v = SparseVec(v)
        return unite!(parent, v, zero(v)) == 0
    end

    inSpan!(x)
    list = [x]
    for y in list, a in aaa
        z = under(y, a)
        inSpan!(z) || push!(list, z)
    end
    return list
end

end # module
