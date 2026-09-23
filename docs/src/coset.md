# Cosets

!!! note
    This is an on-demand module. Load it with `using OrbitAl.coset`.

This module provides a `Coset` type and orbit-based enumeration of right cosets
of a subgroup. Two cosets are identified by testing whether their representatives
differ by an element of the subgroup.

---

## Types

```@docs
OrbitAl.coset.Coset
```

---

## Enumeration

```@docs
OrbitAl.coset.cosets
```
