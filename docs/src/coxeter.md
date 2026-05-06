# Coxeter Groups

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

## Shapes and Conjugacy

```@docs
OrbitAl.coxeter.shape
OrbitAl.coxeter.shapes
OrbitAl.coxeter.coxeterConjugacyClasses
```
