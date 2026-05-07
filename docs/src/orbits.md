# Orbits

This module provides orbit computation algorithms and standard group actions, with options to record word paths, Schreier trees, coset representatives, stabilizers, and more.

---

## Basic Orbit Computation

```@docs
OrbitAl.orbits.orbit
OrbitAl.orbits.orbitl
```

---

## Actions

Action functions passed as the `under` argument to orbit functions:

```@docs
OrbitAl.orbits.onPoints
OrbitAl.orbits.onRight
OrbitAl.orbits.onPairs
OrbitAl.orbits.onSets
OrbitAl.orbits.onWords
```

---

## Orbit Variants

These return richer data structures tracking more information about the orbit:

```@docs
OrbitAl.orbits.orbit_with_words
OrbitAl.orbits.orbit_with_dist
OrbitAl.orbits.orbit_with_tree
OrbitAl.orbits.orbit_with_transversal
OrbitAl.orbits.orbit_with_stabilizer
OrbitAl.orbits.orbit_with_edges
OrbitAl.orbits.orbit_with_images
```

---

## Multiple Starting Points

For cases where orbit computation begins from a set of starting points:

```@docs
OrbitAl.orbits.orbitx
OrbitAl.orbits.orbitx_with_words
OrbitAl.orbits.orbitx_with_dist
OrbitAl.orbits.orbitx_with_edges
```

---

## Utility

```@docs
OrbitAl.orbits.edges_from_images
```

---

## Data Type: `Orbit`

```@docs
OrbitAl.orbits.Orbit
```
