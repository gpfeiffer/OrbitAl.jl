# OrbitAl.jl

Welcome to the documentation for `OrbitAl.jl`, a Julia package for working with permutations and orbits in a simple and composable way.

## Core modules

Loaded automatically with `using OrbitAl`:

- [Permutations](permutation.md)
- [Orbits](orbits.md)
- [Permutation Groups](permgroup.md)
- [Coxeter Groups](coxeter.md)

## On-demand modules

Load explicitly, e.g. `using OrbitAl.syt`:

- [Schreier-Sims Groups](simsgroup.md) — `using OrbitAl.simsgroup`
- [Standard Young Tableaux](syt.md) — `using OrbitAl.syt`
- [Cosets](coset.md) — `using OrbitAl.coset`
- [Involutions](involution.md) — `using OrbitAl.involution`

## Installation

You can install the package from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/gpfeiffer/OrbitAl.jl")
```

## Features

### Core

- **Permutations** — `Perm` type with full arithmetic: composition `*`, inversion `inv`, power `^`, conjugation, cycle decomposition, sign, order.
- **Orbit engine** — BFS-based orbit algorithms in 10+ variants: with words, transversals, stabilizers, edges, images, and multi-seed (`orbitx`) forms.
- **Standard actions** — `onPoints`, `onRight`, `onSets`, `onPairs`, `onWords` ready to use or compose.
- **Permutation groups** — `PermGp` supporting element enumeration, conjugacy classes, subgroup enumeration, membership testing, random element sampling, and intersection.
- **Coxeter groups** — `CoxeterGp` built from a Cartan matrix: root systems, reflections, Coxeter length, reduced words, parabolic subgroups and transversals, conjugacy classes.
- **Visualization** — D3.js force-directed Cayley graphs rendered in Jupyter notebooks.

### On-demand

- **Schreier-Sims groups** — `SimsGp` with cached stabilizer chain for fast repeated membership tests and size computation without enumerating all elements.
- **Standard Young tableaux** — partitions, Newton sums and differences, composition/subset conversions, tableau paths.
- **Cosets** — `Coset` type and `cosets` for orbit enumeration of right cosets of a subgroup.
- **Involutions** — actions and orbit algorithms for involutions and their conjugacy classes in Coxeter groups.

## Usage

### Permutations

```julia
using OrbitAl

p = Perm([2, 3, 1])          # 3-cycle
q = inv(p)
p * q == one(p)              # true
shape(p)                     # [3]
order(p)                     # 3
```

### Orbits

```julia
using OrbitAl

s = Perm([2, 1, 3, 4])       # transposition (1 2)
t = Perm([2, 3, 4, 1])       # 4-cycle (1 2 3 4)

orbit([s, t], 1, onPoints)   # [1, 2, 3, 4]
```

### Coxeter groups

```julia
using OrbitAl

W = CoxeterGp(cartanMat("A", 3))   # symmetric group S4 as a Coxeter group
sizeOfGroup(W)                      # 24
length(coxeterConjugacyClasses(W))  # 5 (partitions of 4)
```

### On-demand: Schreier-Sims

```julia
using OrbitAl.simsgroup

size(cube)   # order of Rubik's cube group
```
