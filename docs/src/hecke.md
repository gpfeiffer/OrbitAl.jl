# Hecke Algebras

!!! note
    This is an on-demand module. Load it with `using OrbitAl.hecke`.

The Iwahori-Hecke algebra ``H`` of a finite Coxeter group ``W`` with
parameter ``q`` has a basis ``T_w``, ``w \in W``, with
``T_w T_s = T_{ws}`` if ``\ell(ws) > \ell(w)``, and the quadratic relations
``T_s^2 = (q-1) T_s + q``.  Here, ``q`` is an element of a field, such as
`2//1`, and ``H`` is a finite-dimensional algebra over that field.  A
`HeckeAlg` does not store the elements of ``W``: products are computed from
reduced words and the action of ``W`` on its roots, so that ``W`` can be as
large as ``E_8``.

An element [`HElt`](@ref OrbitAl.hecke.HElt) of ``H`` knows its algebra.  So
`zero(h)` and `one(h)` work for an element `h`, and `zero(H)` and `one(H)` for
the algebra `H`, but there is no `zero(HElt)`: the type alone does not say
which algebra it means.  Elements support `+`, `-`, `*`, multiplication by
scalars, `==`, and, for the monomial units ``c T_w``, `inv`.

Hecke algebra elements can be coefficients of [sparse
vectors](@ref Sparse-Vectors), acting from the left, and of linear
[Union-Find](@ref).  There, the method of
[`isunit`](@ref OrbitAl.hecke.isunit(::OrbitAl.hecke.HElt)) for `HElt`s
makes `unite!` solve relations for monomial unit coefficients only.  Other
units, such as ``1 - T_s`` for ``q \neq 1``, are not recognised, and a
relation with no monomial unit coefficient waits.

```jldoctest
julia> using OrbitAl.coxeter, OrbitAl.hecke

julia> H = HeckeAlg(CoxeterGp(cartanMat("A", 2)), 2//1)
HeckeAlg(rank 2, q = 2//1)

julia> T1, T2 = Tw(H, 1), Tw(H, 2);

julia> T1 * T1
2 + T1

julia> T1 * T2 * T1 == Tw(H, 1, 2, 1)
true

julia> inv(T1)
-1//2 + 1//2T1

julia> isunit(T1 + one(H))
false
```

---

## Algebras and Elements

```@docs
OrbitAl.hecke.HeckeAlg
OrbitAl.hecke.HElt
OrbitAl.hecke.Tw
OrbitAl.hecke.isunit(::OrbitAl.hecke.HElt)
OrbitAl.hecke.inv(::OrbitAl.hecke.HElt)
```
