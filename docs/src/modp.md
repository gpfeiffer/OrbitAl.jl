# Modular Arithmetic

!!! note
    This is an on-demand module. Load it with `using OrbitAl.modp`.

`Zp{p}` is the field ``\mathbb{F}_p`` of residues modulo a prime `p`, with
`+`, `-`, `*`, `/`, `^` and `inv` (by Fermat's little theorem).  The modulus
is a type parameter, so a vector or matrix over `Zp{p}` carries its modulus
in its type, and integers are promoted to `Zp{p}` where needed.  Products go
through `Int128`, so `p` can be as large as ``2^{62}``.

Since `Zp{p}` is a `Number` with `zero`, `one` and `inv`, it can be the
coefficient type of [sparse vectors](@ref Sparse-Vectors), of linear
[Union-Find](@ref), and of [spinning](@ref Linear-Orbit-Algorithms), without
any changes there.

```jldoctest
julia> using OrbitAl.modp

julia> a = Zp{7}(3)
3 mod 7

julia> inv(a), a^6, a + 5
(5 mod 7, 1 mod 7, 1 mod 7)

julia> using OrbitAl, OrbitAl.linear

julia> A = Zp{7}.([0 1 0; 0 0 1; 1 0 0]);

julia> length(spinning([A], Zp{7}.([1 1 1]), onRight))
1
```

---

## Type and Functions

```@docs
OrbitAl.modp.Zp
```
