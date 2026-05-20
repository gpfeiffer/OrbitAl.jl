# Schreier-Sims Groups

This module provides an alternative permutation group representation built on the
Schreier-Sims stabilizer chain. Unlike `PermGp`, which recomputes the orbit-stabilizer
decomposition on every query, `SimsGp` caches the chain on first use, giving fast
repeated membership tests and size computations without enumerating all elements.

---

## Types

```@docs
OrbitAl.simsgroup.SimsGp
```

---

## Orbit and Stabilizer

```@docs
OrbitAl.simsgroup.orbit_sims
```

---

## Size, Membership and Random Elements

`SimsGp` extends three functions from `Base`:

- **`size(G)`** — returns the group order as a `BigInt`, computed as the product of
  orbit lengths along the stabilizer chain.
- **`in(g, G)`** — tests whether permutation `g` belongs to `G` by sifting through
  the Schreier-Sims chain; each level reduces the element by a coset representative.
- **`rand(G)`** — returns a uniformly random element by multiplying a random coset
  representative at each level of the chain.

---

## Intersection

```@docs
Base.intersect(::OrbitAl.simsgroup.SimsGp, ::OrbitAl.simsgroup.SimsGp)
```
