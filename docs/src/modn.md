# Modular Arithmetic

!!! note
    This is an on-demand module. Load it with `using OrbitAl.modn`.

`Zn{n}` holds the residues modulo `n`, with `+`, `-`, `*`, `/`, `^` and
`inv`.  For a prime `n`, this is the field ``\mathbb{F}_n``.  Otherwise it is
the ring ``\mathbb{Z}/n\mathbb{Z}``: `inv` (by the extended Euclidean algorithm)
raises a `DomainError` for a residue that is not a unit, and `isunit` tells
which residues are units, so that linear Union-Find solves relations for unit
coefficients only.

The modulus is a type parameter, so a vector or matrix over `Zn{n}` carries
its modulus in its type, and integers are promoted to `Zn{n}` where needed.
Products go through `Int128`, so `n` can be as large as ``2^{62}``.

Since `Zn{n}` is a `Number` with `zero`, `one` and `inv`, it can be the
coefficient type of [sparse vectors](@ref Sparse-Vectors), of linear
[Union-Find](@ref), and of [spinning](@ref Linear-Orbit-Algorithms), without
any changes there.

```jldoctest
julia> using OrbitAl.modn

julia> a = Zn{7}(3)
3 mod 7

julia> inv(a), a^6, a + 5
(5 mod 7, 1 mod 7, 1 mod 7)

julia> using OrbitAl, OrbitAl.linear

julia> A = Zn{7}.([0 1 0; 0 0 1; 1 0 0]);

julia> length(spinning([A], Zn{7}.([1 1 1]), onRight))
1
```

---

## Type and Functions

```@docs
OrbitAl.modn.Zn
OrbitAl.modn.isunit(::OrbitAl.modn.Zn)
```
