# Prefix Algebras

!!! note
    This is an on-demand module. Load it with `using OrbitAl.palg`.

A [`PAlg`](@ref OrbitAl.palg.PAlg) is an algebra ``A``, free of rank ``d``
over a commutative ring ``\Lambda``, given by generators ``p_1, \dots, p_r``
and a basis of positive words in them,

```math
b_1 = 1, \quad \dots, \quad b_d, \qquad b_x = b_{\mathrm{parent}(x)} T_{\mathrm{letter}(x)},
```

closed under prefixes. What is stored is the sparse multiplication by one
generator on each side, ``b_x T_j`` and ``T_j b_x``; every other product is
built from those. The algebra need not commute, so the two tables are
different data.

Each generator has a minimal polynomial, `params[j] = [c_0, …, c_{k-1}]` with
the **constant term first**,

```math
T_j^k = c_{k-1} T_j^{k-1} + \dots + c_1 T_j + c_0,
```

and ``c_0`` is a unit, which is what makes ``T_j^{-1}`` an element of ``A``.
Generators that are to share a minimal polynomial share a vector.

The coefficient type carries the ground ring — `Zn{p}` for a run over
``\mathbb{F}_p``, a polynomial type for a symbolic one — and nothing in the
module branches on which. Zero and one are taken from a parameter instance,
not from `zero(T)`, so a coefficient type that defines them only on instances
still works.

## Against `HeckeAlg`

Both are algebras with a basis indexed by words, and they serve opposite
sizes. [`HeckeAlg`](@ref Hecke-Algebras) stores no elements and computes
products from reduced words and the root action, so ``W`` can be as large as
``E_8``; it is the Iwahori-Hecke algebra of a Coxeter group and the quadratic
relation is the same for every generator. A `PAlg` stores both multiplication
tables, so its basis has to be held in memory — but the generators may have
minimal polynomials of different degrees, the group behind it need not be a
Coxeter group, and a product is a table lookup rather than a walk.

## The product

`mul(A, x, y)` expands along the prefixes of the words in the support of `y`,
or along the suffixes of the words in the support of `x`, whichever is
cheaper. The memo then costs one generator multiplication per word reached
rather than one per pair of basis elements, which is what makes a rank in the
tens of thousands workable.

Products accumulate on a dense buffer and become an element once, at the end,
so the cost is the support and not the rank. The buffers are **pooled and
indexed by depth**, because a product may be computed while a caller holds
one — `mul` itself does that along its memo. A single shared buffer would be
silently corrupted and the answer would look plausible. Nothing here is
threaded and the pool does not pretend to be thread-safe.

## The cyclic case

`cyclicAlg(params)` is the one-generator case: the basis is
``1, u, \dots, u^{d-1}``, the algebra is commutative, the two multiplication
tables coincide, and the product is a convolution followed by a fold through
the minimal polynomial. [`ugen`](@ref OrbitAl.palg.ugen),
[`uinv`](@ref OrbitAl.palg.uinv), [`upow`](@ref OrbitAl.palg.upow) and
[`squareParams`](@ref OrbitAl.palg.squareParams) belong to it alone and
assert [`iscyclic`](@ref OrbitAl.palg.iscyclic), so a general algebra handed
to a caller that wants ``u^k`` is an error rather than a wrong answer.

## Units

A basis element is a product of generators, hence a unit, and so is a product
of one with the inverse of another.
[`findUnitInverse`](@ref OrbitAl.palg.findUnitInverse) recognises ``c b``,
``c b^{-1} b'`` and ``c b b'^{-1}``, and verifies every candidate by
multiplying it out. It is the test a search for an invertible pivot wants
where the norm through `multMatrix` is correct and too expensive — the
determinant of a ``d \times d`` matrix, at the ranks this module is for.

Which coefficients count as units of ``\Lambda`` depends on the ground ring,
so the test is an argument: `!iszero` over a field,
[`isunit`](@ref OrbitAl.modn.isunit) over `Zn{n}`, something narrower where
only part of ``\Lambda`` is inverted.

```jldoctest
julia> using OrbitAl.palg

julia> A = cyclicAlg([3//1, 5//1])        # u^2 = 5u + 3
Z2

julia> mul(A, ugen(A), ugen(A))
(3)b1 + (5)b2

julia> mul(A, uinv(A), ugen(A)) == one(A)
true

julia> tau(A, one(A))
1//1

julia> squareParams([3//1, 5//1])         # the minimal polynomial of u^2
2-element Vector{Rational{Int64}}:
 -9
 31
```

---

## Algebras and Elements

```@docs
OrbitAl.palg.PAlg
OrbitAl.palg.PElt
OrbitAl.palg.PAlg(::AbstractString, ::Vector{Vector{T}}, ::Vector{Vector{Int}}, ::Any, ::Any) where T
OrbitAl.palg.cyclicAlg
OrbitAl.palg.iscyclic
```

## Elements and products

```@docs
OrbitAl.palg.basis
OrbitAl.palg.element
OrbitAl.palg.basisInverse
OrbitAl.palg.mul
OrbitAl.palg.rmulgen
OrbitAl.palg.lmulgen
OrbitAl.palg.rmulword
OrbitAl.palg.lmulword
OrbitAl.palg.mulinv
OrbitAl.palg.addall
OrbitAl.palg.sumall
OrbitAl.palg.withAcc
```

## Coordinates and the trace

```@docs
OrbitAl.palg.tau
OrbitAl.palg.coefficient
OrbitAl.palg.coeffs
OrbitAl.palg.fromcoeffs
OrbitAl.palg.multMatrix
```

## Units and the cyclic case

```@docs
OrbitAl.palg.findUnitInverse
OrbitAl.palg.monomialOf
OrbitAl.palg.ugen
OrbitAl.palg.uinv
OrbitAl.palg.upow
OrbitAl.palg.squareParams
```
