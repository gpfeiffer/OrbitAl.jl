# Schreier-Sims Groups

!!! note
    This is an on-demand module. Load it with `using OrbitAl.simsgroup`.

This module is the code of Part 5 of the notebooks: the Schreier-Sims
algorithm, which computes a stabilizer chain with a strong generating set, and
backtrack search through the elements of a group along such a chain.  The
chain is recursive: each `Link` holds the next one, for the stabilizer of its
base point, as in the recursion of `sizeOfGroup` for a `PermGp`.  A `SimsGp`
is a permutation group together with its stabilizer chain, for fast size
computations and membership tests without enumerating its elements.

---

## Stabilizer Chains

```@docs
OrbitAl.simsgroup.Link
OrbitAl.simsgroup.sift
OrbitAl.simsgroup.add!
OrbitAl.simsgroup.schreier_sims
OrbitAl.simsgroup.base
OrbitAl.simsgroup.strong_gens
```

---

## Groups

```@docs
OrbitAl.simsgroup.ASimsGp
OrbitAl.simsgroup.SimsGp
OrbitAl.simsgroup.stabChain
```

The functions `sizeOfGroup`, `memberOfGroup` and `randomGroupElement` work on
a stabilizer chain.  For an `ASimsGp`, such as a `SimsGp` or a `CoxeterGp`, they
use its chain, and `size(G)`, `in(g, G)` and `rand(G)` call them.

```@docs
OrbitAl.permgroup.sizeOfGroup(::Nothing)
OrbitAl.permgroup.memberOfGroup(::Union{Nothing, OrbitAl.simsgroup.Link}, ::Any)
OrbitAl.permgroup.randomGroupElement(::Nothing, ::Any)
```

---

## Backtrack Search

```@docs
OrbitAl.simsgroup.backtrack
OrbitAl.simsgroup.subgp_gens
Base.intersect(::OrbitAl.simsgroup.ASimsGp, ::OrbitAl.simsgroup.ASimsGp)
```
