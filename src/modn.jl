#############################################################################
##
#A  modn.jl                                                           OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  The residues modulo n: the coefficient field F_n, if n is prime
##
##  Zn{n} holds a residue in [0, n).  For a prime n, this is the field F_n;
##  otherwise the ring Z/nZ, whose units are the residues coprime to n.
##  Products go through Int128, so n up to 2^62 is safe; n is a type
##  parameter, so a ring over Zn{n} carries its modulus in its type and no
##  operation has to be told the mode.
##
module modn

import ..unionfind: isunit

export Zn, modulus, isunit

"""
    Zn{n}(x)

The residue of the integer `x` modulo `n`.  For a prime `n`, `Zn{n}` is the
field ``\\mathbb{F}_n``; otherwise it is the ring ``\\mathbb{Z}/n\\mathbb{Z}``, where `inv`
and `isunit` know which residues are units.
"""
struct Zn{n} <: Number
    val::Int64

    function Zn{n}(x::Integer) where n
        n isa Int64 || throw(ArgumentError("the modulus $n is not an Int64"))
        1 < n <= 1 << 62 || throw(ArgumentError("the modulus $n is out of range"))
        new{n}(Int64(mod(x, n)))
    end
end

modulus(::Type{Zn{n}}) where n = n
modulus(::Zn{n}) where n = n

Base.zero(::Type{Zn{n}}) where n = Zn{n}(0)
Base.one(::Type{Zn{n}}) where n = Zn{n}(1)
Base.zero(x::Zn) = zero(typeof(x))
Base.one(x::Zn) = one(typeof(x))
Base.iszero(x::Zn) = x.val == 0
Base.isone(x::Zn) = x.val == 1

Base.:(==)(x::Zn{n}, y::Zn{n}) where n = x.val == y.val
Base.hash(x::Zn{n}, h::UInt) where n = hash(x.val, hash(n, h))

Base.:+(x::Zn{n}, y::Zn{n}) where n = Zn{n}(x.val + y.val)
Base.:-(x::Zn{n}) where n = Zn{n}(-x.val)
Base.:-(x::Zn{n}, y::Zn{n}) where n = Zn{n}(x.val - y.val)
Base.:*(x::Zn{n}, y::Zn{n}) where n = Zn{n}(mod(widemul(x.val, y.val), n))

# by the extended Euclidean algorithm: a DomainError if x is no unit
function Base.inv(x::Zn{n}) where n
    iszero(x) && throw(DivideError())
    Zn{n}(invmod(x.val, n))
end

"""
    isunit(x::Zn)

Whether `x` is a unit: whether its residue is coprime to the modulus.  For a
prime modulus, this means `x` is not zero.
"""
isunit(x::Zn{n}) where n = gcd(x.val, n) == 1

Base.:/(x::Zn{n}, y::Zn{n}) where n = x * inv(y)

function Base.:^(x::Zn{n}, k::Integer) where n
    k < 0 && return inv(x)^(-k)
    Zn{n}(powermod(x.val, k, n))
end

Base.convert(::Type{Zn{n}}, x::Integer) where n = Zn{n}(x)
Base.convert(::Type{Zn{n}}, x::Zn{n}) where n = x
Base.promote_rule(::Type{Zn{n}}, ::Type{<:Integer}) where n = Zn{n}
Zn{n}(x::Zn{n}) where n = x

Base.show(io::IO, x::Zn{n}) where n = print(io, x.val, " mod ", n)

end # module
