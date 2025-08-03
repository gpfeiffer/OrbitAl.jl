# OrbitAl.jl

Welcome to the documentation for `OrbitAl.jl`, a Julia package for working with permutations and orbits in a simple and composable way.

This documentation provides an overview of its modules and functionality.

- [Permutations](permutation.md)
- [Permutation Groups](permgroup.md)
- [Orbits](orbits.md)
- [Coxeter Groups](coxeter.md)

## Installation

You can install the package from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/gpfeiffer/OrbitAl.jl")
```

## Usage
```julia
using OrbitAl

p = Perm([2, 3, 1])
p^3 == one(p)  # true
```

## Features

