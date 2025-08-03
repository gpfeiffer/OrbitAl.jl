# Permutations

This module defines the core `Perm` type and basic operations on permutations in `OrbitAl.jl`.

---

## 📦 Data Type: `Perm`

```@docs
# OrbitAl.permutation.Perm
```

---

## 🏗 Constructors

```@docs
# OrbitAl.permutation.Perm
# OrbitAl.permutation.one
# OrbitAl.permutation.transposition
# OrbitAl.permutation.transpositions
# OrbitAl.permutation.shuffle!
```

---

## 📏 Attributes and Properties

```@docs
# OrbitAl.permutation.degree
# OrbitAl.permutation.domain
# OrbitAl.permutation.shape
# OrbitAl.permutation.order
# OrbitAl.permutation.sign
# OrbitAl.permutation.last_moved
# OrbitAl.permutation.isidentity
```

---

## 🔄 Group Operations

```@docs
# Base.==
# Base.hash
# Base.isless
# Base.inv
# Base.*
# Base./
# Base.^
```

---

## 🎯 Actions

```@docs
# OrbitAl.permutation.^
# OrbitAl.permutation.permuted
```

---

## 🔁 Cycles

```@docs
# OrbitAl.permutation.cycles
```
