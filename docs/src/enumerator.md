# Coset Enumeration

!!! note
    This is an on-demand module. Load it with `using OrbitAl.enumerator`.
    Example presentations are in `OrbitAl.presentations`.

A simple Todd-Coxeter coset enumerator, put together from three ingredients:

- the **orbit algorithm**: each coset is finalized under each generator, while
  new cosets are being defined;
- **Union-Find** (see [`OrbitAl.unionfind`](@ref Union-Find)): cosets that turn
  out to be equal are merged in a forest on the cosets.  Merging two cosets
  also merges their rows, which can cause further coincidences: Union-Find
  with consequences;
- the **variants of the relations**: to find the image of a coset under a
  generator, the relations are rewritten so as to express the generator as a
  word in the others.

A presentation is a named tuple with fields `gens` (the generators `1:k`,
closed under inverses), `invr` (the inverse of each generator), `rels` (the
relations, as pairs of words) and, optionally, `sbgp` (words generating a
subgroup).

```jldoctest
julia> using OrbitAl, OrbitAl.presentations, OrbitAl.enumerator

julia> T = coset_table(presentations.G12, presentations.G12.sbgp);

julia> T.active
24

julia> perms(coset_table(presentations.A2, presentations.A2.sbgp))
2-element Vector{Perm}:
 Perm([2, 1, 3])
 Perm([1, 3, 2])

julia> sizeOfGroup(PermGp(coset_table(presentations.A3, Vector{Int}[])))
24
```

---

## Coset Tables

```@docs
OrbitAl.enumerator.CosetTable
OrbitAl.enumerator.coset_table
OrbitAl.enumerator.is_active(::OrbitAl.enumerator.CosetTable, ::Any)
OrbitAl.enumerator.active_cosets(::OrbitAl.enumerator.CosetTable)
OrbitAl.enumerator.perms
```
