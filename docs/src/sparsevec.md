# Sparse Vectors

!!! note
    This is an on-demand module. Load it with `using OrbitAl.sparsevec`.

A `SparseVec` stores the positions `poss` of the nonzero coefficients of a
vector, in increasing order, together with the coefficients `vals`.  It has no
length, so vectors can involve more and more basis vectors, such as the cosets
of an enumeration.  The last position, which linear Union-Find asks for all the
time, is simply `v.poss[end]`.

Sparse vectors support `+`, `-`, multiplication by scalars, division by
invertible scalars, indexing `v[i]`, `zero`, `iszero`, `length` (the number of
nonzero coefficients), `eltype` (the type of the coefficients), `==` and
`hash`.

```jldoctest
julia> using OrbitAl.sparsevec

julia> v = SparseVec(3 => 2//1, 1 => 1//1)
e1 + (2)e3

julia> v[3], v[2]
(2//1, 0//1)

julia> v - 2 * unitVec(Rational{Int}, 3)
e1

julia> SparseVec(Rational{Int}[0 1 0 4])
e2 + (4)e4
```

---

## Type and Constructors

```@docs
OrbitAl.sparsevec.SparseVec
OrbitAl.sparsevec.unitVec
```
