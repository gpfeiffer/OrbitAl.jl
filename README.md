# `OrbitAl.jl`

A Julia package for working with permutations and orbits in a simple and composable way.

## Features

- **Permutations** — `Perm` type with full arithmetic: composition `*`, inversion `inv`, power `^`, conjugation, cycle decomposition, sign, order.
- **Orbit engine** — BFS-based orbit algorithms in 10+ variants: with words, transversals, stabilizers, edges, images, and multi-seed (`orbitx`) forms.
- **Standard actions** — `onPoints`, `onRight`, `onSets`, `onPairs`, `onWords` ready to use or compose.
- **Permutation groups** — `PermGp` supporting element enumeration, conjugacy classes, subgroup enumeration, membership testing, and random element sampling.
- **Coxeter groups** — `CoxeterGp` built from a Cartan matrix: root systems, reflections, Coxeter length, reduced words, parabolic subgroups and transversals, conjugacy classes.
- **Standard Young tableaux** — partitions, Newton polynomials, composition/subset conversions, tableau paths.
- **Visualization** — D3.js force-directed Cayley graphs rendered in Jupyter notebooks.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/gpfeiffer/OrbitAl.jl")
```

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

## Documentation

Full documentation is available at <https://gpfeiffer.github.io/OrbitAl.jl/dev/>.
