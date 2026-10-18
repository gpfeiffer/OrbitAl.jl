# Coxeter Groups

!!! note
    This is an on-demand module. Load it with `using OrbitAl.coxeter`.

This module constructs finite Coxeter groups as permutation groups acting on their
root systems, and provides algorithms for lengths, reduced words, reflections,
parabolic subgroups, and conjugacy classes.

---

## Building a Coxeter Group

```@docs
OrbitAl.coxeter.coxeterGraph
OrbitAl.coxeter.coxeterMat
OrbitAl.coxeter.cartanMat
OrbitAl.coxeter.CoxeterGp
```

---

## Length and Words

```@docs
OrbitAl.coxeter.coxeterLength
OrbitAl.coxeter.coxeterWord
OrbitAl.coxeter.permCoxeterWord
OrbitAl.coxeter.reflections
```

---

## Parabolic Subgroups

```@docs
OrbitAl.coxeter.longestElt
OrbitAl.coxeter.parabolicTransversal
OrbitAl.coxeter.prefixes
OrbitAl.coxeter.prefixes_with_edges
```

---

## Poincaré Polynomials

Polynomials are lists of coefficients, constant term first.

```@docs
OrbitAl.coxeter.qint
OrbitAl.coxeter.conv
OrbitAl.coxeter.lengthPoly
OrbitAl.coxeter.poincareByChain
OrbitAl.coxeter.coxeterNumber
```

---

## Shapes and Conjugacy

```@docs
OrbitAl.coxeter.shape
OrbitAl.coxeter.shapes
OrbitAl.coxeter.coxeterConjugacyClasses
```
