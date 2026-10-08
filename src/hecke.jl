#############################################################################
##
#A  hecke.jl                                                          OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  Iwahori-Hecke algebras of finite Coxeter groups, in the basis T_w
##
module hecke

using ..permutation
using ..coxeter
import ..coxeter: isLeftDescent
import ..unionfind: isunit

export HeckeAlg, HElt, Tw, isunit

"""
    HeckeAlg(W, q)

The Iwahori-Hecke algebra of the finite Coxeter group `W`, with parameter
`q`, an element of a field, such as `2//1`.  It has a basis ``T_w``,
``w \\in W``, with ``T_w T_s = T_{ws}`` if ``\\ell(ws) > \\ell(w)``, and
``T_s^2 = (q-1) T_s + q``.

A Hecke algebra makes its own elements: `zero(H)`, `one(H)` and
[`Tw`](@ref)`(H, s...)`.  It does not store the elements of `W`, so `W` can
be large.
"""
struct HeckeAlg{C}
    W::CoxeterGp
    q::C
end

Base.show(io::IO, H::HeckeAlg) =
    print(io, "HeckeAlg(rank ", length(H.W.gens), ", q = ", H.q, ")")

"""
    HElt

An element of a Hecke algebra `H`: the nonzero coefficients of the basis
elements ``T_w``, as a dictionary from `w` to its coefficient.  Each element
knows its algebra `H`, so `zero(h)` and `one(h)` work, but `zero(HElt)` does
not.
"""
struct HElt{C}
    H::HeckeAlg{C}
    coeffs::Dict{Perm, C}      # the coefficient of T_w, nonzero only
end

##  drop the zero coefficients
clean(H, d) = HElt(H, filter(p -> !iszero(p.second), d))

"""
    Tw(H, s...)

The basis element ``T_w`` of `H`, for the word `s...` of `w` in the
generators of ``W``.  `Tw(H)` is the identity.
"""
Tw(H::HeckeAlg, word::Int...) = HElt(H, Dict(permCoxeterWord(H.W, [word...]) => one(H.q)))

Base.zero(H::HeckeAlg{C}) where C = HElt(H, Dict{Perm, C}())
Base.one(H::HeckeAlg) = Tw(H)
Base.zero(h::HElt) = zero(h.H)
Base.one(h::HElt) = one(h.H)
Base.iszero(h::HElt) = isempty(h.coeffs)
Base.isone(h::HElt) = h == one(h)

Base.:(==)(h::HElt, g::HElt) = h.H === g.H && h.coeffs == g.coeffs
Base.hash(h::HElt, x::UInt) = hash(h.coeffs, x)
Base.broadcastable(h::HElt) = Ref(h)

##  the algebra of h and g, which have to be the same
function algebra(h::HElt, g::HElt)
    h.H === g.H || error("elements of different Hecke algebras")
    return h.H
end

Base.:+(h::HElt, g::HElt) = clean(algebra(h, g), mergewith(+, h.coeffs, g.coeffs))
Base.:*(c::Number, h::HElt) = clean(h.H, Dict(w => c * a for (w, a) in h.coeffs))
Base.:-(h::HElt) = (-1) * h
Base.:-(h::HElt, g::HElt) = h + (-g)

##  whether l(ws) < l(w), from the action of w on the roots
isRightDescent(W, w, s) = isLeftDescent(W, inv(w), s)

##  h * T_s
function mulgen(h::HElt, s::Int)
    H = h.H
    d = empty(h.coeffs)
    for (w, c) in h.coeffs
        ws = w * H.W.gens[s]
        if !isRightDescent(H.W, w, s)
            d[ws] = get(d, ws, 0) + c
        else
            d[w] = get(d, w, 0) + (H.q - 1) * c
            d[ws] = get(d, ws, 0) + H.q * c
        end
    end
    return clean(H, d)
end

function Base.:*(h::HElt, g::HElt)
    H = algebra(h, g)
    terms = (c * foldl(mulgen, coxeterWord(H.W, w); init = h) for (w, c) in g.coeffs)
    return sum(terms; init = zero(H))
end

##  h * T_s⁻¹, where T_s⁻¹ = q⁻¹ T_s + (q⁻¹ - 1)
mulinv(h::HElt, s::Int) = (p = inv(h.H.q); p * mulgen(h, s) + (p - 1) * h)

"""
    isunit(h::HElt)

Whether `h` is a monomial unit ``c T_w``: a single term, with ``c`` a unit,
and ``w = 1`` or ``q`` a unit.

Units that are not monomials, such as ``1 - T_s`` for ``q \\neq 1``, are not
recognised.  Linear [Union-Find](@ref) then lets a relation with such a
coefficient wait, rather than solve it.
"""
function isunit(h::HElt)
    length(h.coeffs) == 1 || return false
    w, c = only(h.coeffs)
    return isunit(c) && (isidentity(w) || isunit(h.H.q))
end

"""
    inv(h::HElt)

The inverse of the monomial unit ``h = c T_w``: if ``w = s_1 \\cdots s_k`` is
reduced, then ``h^{-1} = c^{-1} T_{s_k}^{-1} \\cdots T_{s_1}^{-1}``, where
``T_s^{-1} = q^{-1} T_s + (q^{-1} - 1)``.
"""
function Base.inv(h::HElt)
    isunit(h) || error("not a monomial unit: ", h)
    w, c = only(h.coeffs)
    return foldl(mulinv, reverse(coxeterWord(h.H.W, w)); init = inv(c) * one(h))
end

function Base.show(io::IO, h::HElt)
    iszero(h) && return print(io, "0")
    terms = sort([(coxeterWord(h.H.W, w), c) for (w, c) in h.coeffs], by = t -> (length(t[1]), t[1]))
    cf(c) = sprint(show, c; context = :typeinfo => typeof(c))
    term(word, c) = isempty(word) ? cf(c) : (isone(c) ? "" : c == -1 ? "-" : cf(c)) * "T" * join(word)
    print(io, join([term(word, c) for (word, c) in terms], " + "))
end

end # module
