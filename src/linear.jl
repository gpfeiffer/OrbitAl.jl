#############################################################################
##
#A  linear.jl                                                         OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  Linear orbit algorithms: spinning, and spinning with images
##
module linear

using ..sparsevec
using ..unionfind
import ..unionfind: lastDead

export spinning, spinning_with_images

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

"""
    spinning_with_images(aaa, x, under)

Spinning, with images: returns the basis `list` that `spinning` returns, and
for each operator `a` in `aaa` the list of the images `under(y, a)` of the
basis vectors `y`, as `SparseVec`s of coordinates in terms of `list`.  These
are the rows of the matrices of the operators on the subspace.

The coordinates come from linear Union-Find with witnesses: alongside each
relation `e_i - parent[i]`, `coeffs[i]` records it as a combination of the
basis vectors.  Reducing a vector to `0` then expresses it in terms of the
basis, as in extended Gaussian elimination.
"""
function spinning_with_images(aaa, x, under)
    T = eltype(x)

    # the span of list, as relations, and the relations as combinations of list
    parent = Dict{Int, SparseVec{T}}()
    coeffs = Dict{Int, SparseVec{T}}()

    # the coordinates of z in terms of list, after adding z to list if it is new
    function coordinates!(z)
        v, c = SparseVec(z), zero(SparseVec{T})
        while (k = lastDead(parent, v)) > 0        # find, collecting coefficients
            i, a = v.poss[k], v.vals[k]
            v -= a * (unitVec(T, i) - parent[i])
            c += a * coeffs[i]
        end
        length(v) > 0 || return c                  # z = sum of c[l] list[l]
        push!(list, z)
        m = length(list)
        i, g = v.poss[end], v.vals[end]            # v = z - sum of c[l] list[l]
        parent[i] = unitVec(T, i) - v / g
        coeffs[i] = (unitVec(T, m) - c) / g
        return unitVec(T, m)
    end

    list = typeof(x)[]
    coordinates!(x)
    images = [SparseVec{T}[] for _ in aaa]
    for y in list, (k, a) in enumerate(aaa)
        push!(images[k], coordinates!(under(y, a)))
    end
    return (list = list, images = images)
end

end # module
