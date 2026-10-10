# Schreier-Sims Groups

!!! note
    This is an on-demand module. Load it with `using OrbitAl.simsgroup`.

This module is the code of Part 5 of the notebooks: the Schreier-Sims
algorithm, which computes a stabilizer chain with a strong generating set, and
backtrack search through the elements of a group along such a chain.  A
`SimsGp` is a permutation group whose stabilizer chain is computed on first
use, for fast size computations and membership tests without enumerating its
elements.

---

## Stabilizer Chains

```@docs
OrbitAl.simsgroup.Link
OrbitAl.simsgroup.sift
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
OrbitAl.permgroup.sizeOfGroup(::Vector{OrbitAl.simsgroup.Link})
OrbitAl.permgroup.memberOfGroup(::Vector{OrbitAl.simsgroup.Link}, ::Any)
OrbitAl.permgroup.randomGroupElement(::Vector{OrbitAl.simsgroup.Link}, ::Any)
```

---

## Backtrack Search

```@docs
OrbitAl.simsgroup.backtrack
OrbitAl.simsgroup.subgp_gens
Base.intersect(::OrbitAl.simsgroup.ASimsGp, ::OrbitAl.simsgroup.ASimsGp)
```
