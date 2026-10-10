# Permutation Groups

This module provides lightweight permutation groups built on top of the orbit algorithms.
A group is defined by a list of generators and an identity element.  The abstract type
`APermGp` asks its subtypes for `sizeOfGroup`, `memberOfGroup` and `randomGroupElement`,
and builds the generic algorithms (elements, conjugacy classes, subgroups, intersections)
on them.  `PermGp` computes these three by orbit computations, as in Part 1 of the
notebooks; `SimsGp` and `CoxeterGp` use a stabilizer chain instead.

---

## Types

```@docs
OrbitAl.permgroup.APermGp
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
