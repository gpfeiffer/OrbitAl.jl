# Standard Young Tableaux

This module provides tools for working with integer partitions and standard
Young tableaux, built on top of the orbit machinery.

---

## Newton Sums and Differences

Mutually inverse transformations on integer vectors:

```@docs
OrbitAl.syt.newtonDif
OrbitAl.syt.newtonSum
OrbitAl.syt.newtonDifR
OrbitAl.syt.newtonSumR
```

---

## Compositions and Subsets

Bijection between compositions of `n` and subsets of `1:n-1`:

```@docs
OrbitAl.syt.compositionSubset
OrbitAl.syt.subsetComposition
```

---

## Enumeration

```@docs
OrbitAl.syt.takeAway
OrbitAl.syt.subsets
OrbitAl.syt.partitions
```

---

## Standard Young Tableaux

```@docs
OrbitAl.syt.standardYTs
OrbitAl.syt.tableau_path
```
