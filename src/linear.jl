#############################################################################
##
#A  linear.jl                                                         OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  Linear orbit algorithms: spinning, and spinning with images,
#C  and Norton's irreducibility test (the MeatAxe)
##
module linear

using ..sparsevec
using ..unionfind
import ..orbits: onRight
import ..modn: modulus

export spinning, spinning_with_images, left_nullspace, norton

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
        v = SparseVec(z)
        v, c = find_with_coeffs(parent, coeffs, v, zero(v))
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

"""
    left_nullspace(A)

A basis of the left null space `{d : d A = 0}` of the matrix `A`, as a list
of row vectors.

The rows of `A` are entered one by one into linear Union-Find, as relations,
and `coeffs` records each relation as a combination of the rows, as in
`spinning_with_images`.  A row that reduces to `0` gives a null vector
`d = e_k - c`; its last nonzero entry is at position `k`, so the null vectors
are linearly independent.
"""
function left_nullspace(A)
    n = size(A, 1)
    T = eltype(A)

    # the rows of A, as relations, and the relations as combinations of the rows
    parent = Dict{Int, SparseVec{T}}()
    coeffs = Dict{Int, SparseVec{T}}()

    nulls = Matrix{T}[]
    for k in 1:n
        v = SparseVec(A[k, :])
        v, c = find_with_coeffs(parent, coeffs, v, zero(v))    # v = (e_k - c) A
        if length(v) == 0
            d = unitVec(T, k) - c                                # d A = 0
            push!(nulls, permutedims([d[i] for i in 1:n]))
        else
            i, g = v.poss[end], v.vals[end]        # a new relation, as in unite!
            parent[i] = unitVec(T, i) - v / g
            coeffs[i] = (unitVec(T, k) - c) / g
        end
    end
    return nulls
end

"""
    norton(gens; tries = 50)

Norton's irreducibility test, the core of the MeatAxe, for the module over a
prime field `Zn{p}` given by the matrices `gens`, acting on row vectors.  A
random element θ of the algebra, with a null space of dimension 1 spanned by
`v`, decides: the module is irreducible if and only if `v` spins to the
whole space, and a vector in the null space of θ^T spins to the whole space
under the transposed matrices.

Returns `(irreducible = false, submodule = U)` or `(irreducible = false,
dual_submodule = W)` with a basis of a proper submodule, `(irreducible = true,)`,
or `(irreducible = missing,)` if no suitable θ was found in `tries` attempts.
"""
function norton(gens; tries = 50)
    n = size(gens[1], 1)
    F = eltype(gens[1])
    field = F.(0:modulus(F)-1)
    for _ in 1:tries
        words = [prod(rand(gens, rand(1:4))) for _ in 1:3]
        θ0 = sum(rand(field) * w for w in words)
        for λ in field
            θ = θ0 - λ * one(θ0)
            N = left_nullspace(θ)
            length(N) == 1 || continue
            U = spinning(gens, N[1], onRight)
            length(U) < n && return (irreducible = false, submodule = U)
            W = spinning(transpose.(gens), left_nullspace(transpose(θ))[1], onRight)
            length(W) < n && return (irreducible = false, dual_submodule = W)
            return (irreducible = true,)
        end
    end
    return (irreducible = missing,)
end

end # module
