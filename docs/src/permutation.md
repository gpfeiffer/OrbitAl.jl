# Permutations

This module defines the core `Perm` type and basic operations on permutations in `OrbitAl.jl`.

---

## 📦 Data Type: `Perm`

```@docs
OrbitAl.permutation.Perm
```

---

## 🏗 Constructors

```@docs
OrbitAl.permutation.one
OrbitAl.permutation.transpositions
OrbitAl.permutation.shuffle!
OrbitAl.permutation.rand
```

---

## 📏 Attributes and Properties

```@docs
OrbitAl.permutation.degree
OrbitAl.permutation.domain
OrbitAl.permutation.shape
OrbitAl.permutation.order
OrbitAl.permutation.sign
OrbitAl.permutation.last_moved
OrbitAl.permutation.isidentity
OrbitAl.permutation.cycles
```

---

## 🔄 Group Operations

```@docs
==(::OrbitAl.permutation.Perm, ::OrbitAl.permutation.Perm)
hash(::OrbitAl.permutation.Perm, ::UInt)
isless(::OrbitAl.permutation.Perm, ::OrbitAl.permutation.Perm)
inv(::OrbitAl.permutation.Perm)
*(::OrbitAl.permutation.Perm, ::OrbitAl.permutation.Perm)
/(::OrbitAl.permutation.Perm, ::OrbitAl.permutation.Perm)
^(::OrbitAl.permutation.Perm, ::Int)
```

---

## 🎯 Actions

```@docs
^(::Int, ::OrbitAl.permutation.Perm)
OrbitAl.permutation.permuted
```
