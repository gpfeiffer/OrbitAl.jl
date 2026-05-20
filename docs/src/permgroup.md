# Permutation Groups

This module provides lightweight permutation groups built on top of the orbit algorithms.
A group is defined by a list of generators and an identity element; all group-theoretic
operations (size, membership, conjugacy classes) are derived from orbit computations.

---

## Types

```@docs
OrbitAl.permgroup.PermGp
```

---

## Element Enumeration

```@docs
OrbitAl.permgroup.elements
OrbitAl.permgroup.subgroups
```

---

## Size and Membership

```@docs
OrbitAl.permgroup.sizeOfGroup
OrbitAl.permgroup.memberOfGroup
OrbitAl.permgroup.randomGroupElement
```

---

## Conjugacy

```@docs
OrbitAl.permgroup.conjClasses
OrbitAl.permgroup.subgpClasses
OrbitAl.permgroup.closure
```

---

## Subgroup Structure

```@docs
OrbitAl.permgroup.isPrimePower
OrbitAl.permgroup.zuppos
```

---

## Intersection

```@docs
Base.intersect(::OrbitAl.permgroup.APermGp, ::OrbitAl.permgroup.APermGp)
```
