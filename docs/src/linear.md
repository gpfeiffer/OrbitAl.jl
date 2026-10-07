# Linear Orbit Algorithms

!!! note
    This is an on-demand module. Load it with `using OrbitAl.linear`.

Spinning is the linear version of the orbit algorithm: starting from a vector
`x`, it finds a basis of the smallest subspace that contains `x` and is
invariant under the given operators.  A vector is new if it is not in the span
of the basis so far, which linear Union-Find (see
[`OrbitAl.unionfind`](@ref Union-Find)) tests incrementally.

`spinning_with_images` also expresses the image of each basis vector under
each operator in terms of the basis: the rows of the matrices of the operators
on the subspace.

```jldoctest
julia> using OrbitAl, OrbitAl.linear

julia> A = Rational{Int}[0 1 0; 0 0 1; 1 0 0];   # a 3-cycle, acting on row vectors

julia> vcat(spinning([A], Rational{Int}[1 -1 0], onRight)...)
2×3 Matrix{Rational{Int64}}:
 1  -1   0
 0   1  -1

julia> spinning_with_images([A], Rational{Int}[1 -1 0], onRight).images[1]
2-element Vector{OrbitAl.sparsevec.SparseVec{Rational{Int64}}}:
 e2
 (-1)e1 + (-1)e2
```

---

## Spinning

```@docs
OrbitAl.linear.spinning
OrbitAl.linear.spinning_with_images
```
