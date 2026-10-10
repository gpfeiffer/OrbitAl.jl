# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Run all tests
julia --project -e 'using Pkg; Pkg.test("OrbitAl")'

# Run a single test file interactively
julia --project test/runtests.jl

# Build documentation locally
julia --project=docs docs/make.jl

# Start Julia REPL with the package loaded
julia --project -e 'using OrbitAl'
```

## Architecture

OrbitAl.jl is a Julia package for computational group theory and orbit algorithms. The architecture layers cleanly from primitives to applications:

**Foundation layer**
- `permutation.jl` — `Perm` struct (stored as `Vector{Int}` of images). All permutation arithmetic: `*`, `inv`, `^` (conjugation), cycles, sign, order. The `^` operator is overloaded for point/set/word actions.
- `orbits.jl` — The core BFS orbit engine with ~10 variants: `orbit`, `orbitl`, `orbit_with_words`, `orbit_with_tree`, `orbit_with_transversal`, `orbit_with_stabilizer`, `orbit_with_images`, `orbit_with_edges`. Standard actions (`onPoints`, `onRight`, `onSets`, `onPairs`, `onWords`) live here.

**Group layer** (depends on permutation + orbits)
- `permgroup.jl` — `APermGp`, the abstract type: subtypes provide `sizeOfGroup`, `memberOfGroup`, `randomGroupElement`, and get `size`/`in`/`rand` and the generic algorithms (elements, closure, conjugacy classes, subgroups). `PermGp` struct: those three by the recursive orbit-stabilizer method of Part 1. Also `isPrimePower`, `zuppos`, and `intersect` (extending `Base.intersect`/`∩`) for orbit-stabilizer group intersection. Built on top of orbit algorithms.
- `simsgroup.jl` — the code of Part 5 of the notebooks: `Link` (a recursive chain: each link holds the `next`, `nothing` at the end), `sift`, `add!`, `schreier_sims` (stabilizer chains with a strong generating set), `backtrack` and `subgp_gens` (backtrack search), `ASimsGp <: APermGp`, the abstract type of groups with a stabilizer chain (`stabChain`), on which `sizeOfGroup`, `memberOfGroup`, `randomGroupElement` and `intersect` use the chain, and `SimsGp`, a group together with its chain. Included before `coxeter.jl`. On-demand: `using OrbitAl.simsgroup`.
- `bfsdfs.jl` — `Node` struct, generic BFS/DFS tree traversal utilities.

**Algebraic structures** (depend on group layer)
- `coxeter.jl` — `CoxeterGp` from a Coxeter matrix or graph; `ACoxeterGp <: ASimsGp`, with the stabilizer chain kept in `W.data[:chain]`. Reflection matrices, root systems, conjugacy classes, lengths, parabolic subgroups.
- `involution.jl` — Actions on involutions; conjugacy class enumeration relative to Coxeter groups.
- `shifts.jl` — Cyclic shift enumeration with edge tracking.
- `syt.jl` — Standard Young tableaux: partitions, Newton polynomials, composition↔subset conversions, tableau paths.

**Coset enumeration** (experimental/in-progress)
- `data/presentations.jl` — module `presentations`: group presentations as named tuples (`gens`, `invr`, `rels`, `sbgp`), e.g. `A2`, `G333`, `M12`, and the Fibonacci group `F27`, a hard case for coset enumeration.
- `data/permgroups.jl` — module `permgroups`: generator lists of the Mathieu groups `m11` … `m24` (from GAP's `MathieuGroup`) and the Rubik's cube group `cube`.
- `variants.jl` — Relation variant generation from presentations.
- `enumerator.jl` — Todd-Coxeter-style modular coset enumerator.
- `coset.jl` — `Coset` struct; incomplete coset action implementation.

**Visualization**
- `plotting.jl` — D3.js force-directed graphs rendered in Jupyter notebooks; colored edge support.

## Key Conventions

- Permutations use **1-based indexing** throughout.
- The caret `^` operator is heavily overloaded: `x ^ g` applies permutation `g` to point/set/word `x`, and `g ^ h` conjugates `g` by `h`.
- Orbit functions follow a consistent BFS pattern — the primary `orbit` returns a `Vector`, while `orbitl` returns an `Orbit` struct containing both the list and a position lookup.
- The coset enumeration files (`presentations.jl`, `variants.jl`, `enumerator.jl`, `coset.jl`) are experimental and not yet part of the exported API.

## Documentation

Inline docstrings use `@doc`-compatible format and are tested via Documenter.jl doctests. When adding exported functions, add a docstring with at least one `jldoctest` block. The docs CI runs `julia --project=docs docs/make.jl` on push to `main`.
