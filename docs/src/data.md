# Example Groups

!!! note
    These are on-demand modules, in the folder `src/data`.  Load them with
    `using OrbitAl.presentations` or `using OrbitAl.permgroups`, or refer to
    their contents as `presentations.A3`, `permgroups.m24`, and so on.

---

## Presentations

`OrbitAl.presentations` holds presentations for coset enumeration
(`OrbitAl.enumerator`) and vector enumeration (`OrbitAl.vectorenum`).  Each is
a named tuple with the generators `gens`, their inverses `invr`, the relations
`rels` as pairs of words, and generators of some subgroups `sbgp`:

- the Coxeter groups `A2` and `A3`;
- the exceptional complex reflection groups `G4`, …, `G13`, `G17`, `G24`, `G29`,
  and the imprimitive group `G333` = ``G(3,3,3)``;
- the sporadic simple groups `M12` and `J1`;
- the Fibonacci group `F27` = ``F(2,7)``, cyclic of order 29, a hard case for
  coset enumeration: over the trivial subgroup, `coset_table` defines 123415
  cosets before it finds the 29.

---

## Permutation Groups

`OrbitAl.permgroups` holds lists of generating permutations, to be turned into
a group as needed:

- the Mathieu groups `m11`, `m12`, `m22`, `m23` and `m24`, with the generators
  of `MathieuGroup(n)` in GAP;
- the Rubik's cube group `cube`, on the 48 movable facets of the cube.

```jldoctest
julia> using OrbitAl, OrbitAl.permgroups, OrbitAl.simsgroup

julia> size(SimsGp(permgroups.m24, Perm(24)))
244823040

julia> size(SimsGp(permgroups.cube, Perm(48)))
43252003274489856000
```
