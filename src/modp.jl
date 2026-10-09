#############################################################################
##
#A  modp.jl                                                           OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  The residues modulo p: the coefficient field F_p, if p is prime
##
##  Zp{p} holds a residue in [0, p).  For a prime p, this is the field F_p;
##  otherwise the ring Z/pZ, whose units are the residues coprime to p.  Products go through Int128, so p up to
##  2^62 is safe; p is a type parameter, so a ring over Zp{p} carries its
##  modulus in its type and no operation has to be told the mode.
##
module modp

import ..unionfind: isunit

export Zp, modulus, isunit

"""
    Zp{p}(x)

The residue of the integer `x` modulo `p`.  For a prime `p`, `Zp{p}` is the
field ``\\mathbb{F}_p``; otherwise it is the ring ``\\mathbb{Z}/p\\mathbb{Z}``, where `inv`
and `isunit` know which residues are units.
"""
struct Zp{p} <: Number
    val::Int64

    function Zp{p}(x::Integer) where p
        p isa Int64 || throw(ArgumentError("the modulus $p is not an Int64"))
        1 < p <= 1 << 62 || throw(ArgumentError("the modulus $p is out of range"))
        new{p}(Int64(mod(x, p)))
    end
end

modulus(::Type{Zp{p}}) where p = p
modulus(::Zp{p}) where p = p

Base.zero(::Type{Zp{p}}) where p = Zp{p}(0)
Base.one(::Type{Zp{p}}) where p = Zp{p}(1)
Base.zero(x::Zp) = zero(typeof(x))
Base.one(x::Zp) = one(typeof(x))
Base.iszero(x::Zp) = x.val == 0
Base.isone(x::Zp) = x.val == 1

Base.:(==)(x::Zp{p}, y::Zp{p}) where p = x.val == y.val
Base.hash(x::Zp{p}, h::UInt) where p = hash(x.val, hash(p, h))

Base.:+(x::Zp{p}, y::Zp{p}) where p = Zp{p}(x.val + y.val)
Base.:-(x::Zp{p}) where p = Zp{p}(-x.val)
Base.:-(x::Zp{p}, y::Zp{p}) where p = Zp{p}(x.val - y.val)
Base.:*(x::Zp{p}, y::Zp{p}) where p = Zp{p}(mod(widemul(x.val, y.val), p))

# by the extended Euclidean algorithm: a DomainError if x is no unit
function Base.inv(x::Zp{p}) where p
    iszero(x) && throw(DivideError())
    Zp{p}(invmod(x.val, p))
end

"""
    isunit(x::Zp)

Whether `x` is a unit: whether its residue is coprime to the modulus.  For a
prime modulus, this means `x` is not zero.
"""
isunit(x::Zp{p}) where p = gcd(x.val, p) == 1

Base.:/(x::Zp{p}, y::Zp{p}) where p = x * inv(y)

function Base.:^(x::Zp{p}, n::Integer) where p
    n < 0 && return inv(x)^(-n)
    Zp{p}(powermod(x.val, n, p))
end

Base.convert(::Type{Zp{p}}, x::Integer) where p = Zp{p}(x)
Base.convert(::Type{Zp{p}}, x::Zp{p}) where p = x
Base.promote_rule(::Type{Zp{p}}, ::Type{<:Integer}) where p = Zp{p}
Zp{p}(x::Zp{p}) where p = x

Base.show(io::IO, x::Zp{p}) where p = print(io, x.val, " mod ", p)

end # module
