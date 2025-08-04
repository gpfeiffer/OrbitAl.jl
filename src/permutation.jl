#############################################################################
##
#A  permutation.jl                                                    OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  Permutations, group operations and some standard actions
##
module permutation

import Base: length, hash, inv, isless, one, sign, rand, ==, *, /, ^

export Perm, degree, domain, cycles, shape, order, isidentity, last_moved
export permuted, transpositions, shuffle!

##  Perm data type
"""
    Perm

A permutation represented by a list of integers.
The entry at index `i` indicates the image of `i` under the permutation.

Use `Perm(n)` for the identity permutation of size `n`, or `Perm(n, cycles)` to construct from cycle notation.
"""
struct Perm
    list::Vector{Int}
end

##  Constructors
"""
    Perm(n::Int)

Construct the identity permutation on `1:n`.
"""
Perm(n::Int) = Perm(collect(1:n))

"""
    one(Perm, n::Int)

Alias for `Perm(n)` — the identity permutation on `n` elements.
"""
one(Perm, n::Int) = Perm(n)

"""
    one(perm::Perm)

Identity permutation with the same degree as `perm`.
"""
one(perm::Perm) = Perm(domain(perm))

"""
    Perm(n::Int, cycles::Array)

Construct a permutation of degree `n` from disjoint cycles.
Each cycle is a vector of indices.
"""
function Perm(n::Int, cycles::Array)
    perm = Perm(n)
    for c in cycles
        perm.list[c] = c[(1:end) .% end .+ 1]
    end
    return perm
end

##  Basic attributes
"""
    degree(perm::Perm)

Returns the degree (length of the domain) of the permutation.
"""
degree(perm::Perm) = length(perm.list)

"""
    domain(perm::Perm)

Returns the domain of the permutation as range `1:n`.
"""
domain(perm::Perm) = 1:degree(perm)

##  Equality, hash, comparison
"""
    ==(p::Perm, q::Perm)

Check two permutations for equality.
"""
==(perm::Perm, other::Perm) = perm.list == other.list

"""
    hash(perm::Perm, h::UInt)

Compute the hash of a permutation.
"""
hash(perm::Perm, h::UInt) = hash(perm.list, h)

"""
    isless(p::Perm, q::Perm)

Lexicographic comparison of two permutations.
"""
isless(perm::Perm, other::Perm) = perm.list < other.list

##  Action (on points)

"""
    x ^ p

Apply the permutation `p` to the point or vector `x`.
"""
^(x::Int, perm::Perm) = perm.list[x]
^(xs::Vector{Int}, perm::Perm) = perm.list[xs]

##  Inverse permutation
"""
    inv(perm::Perm)

Returns the inverse of the permutation.
"""
function inv(perm::Perm)
    other = similar(perm.list)
    other[perm.list] = eachindex(other)
    return Perm(other)
end

##  Group operations
"""
    p * q

Permutation multiplication: `i^(p*q) == (i^p)^q`.
"""
*(perm::Perm, other::Perm) = Perm(perm.list^other)

"""
    p / q

Right division: `p * inv(q)`
"""
/(perm::Perm, other::Perm) = perm * inv(other)

"""
    p ^ q

Conjugation of `p` by `q`: `q⁻¹ * p * q`
"""
^(perm::Perm, other::Perm) = inv(other) * perm * other

"""
    p ^ n

Raise the permutation `p` to the integer power `n`.
"""
function ^(perm::Perm, n::Int)
    n == 1 && return perm
    n == 0 && return one(perm)
    n < 0 && return inv(perm)^(-n)
    q, r = divrem(n, 2)
    return perm^q * perm^q * perm^r
end

##  shape, order, sign
"""
    shape(perm::Perm)

Returns the list of cycle lengths, sorted decreasingly.
"""
shape(perm::Perm) = sort(length.(cycles(perm)), rev=true)

"""
    order(perm::Perm)

Returns the order of the permutation.
"""
order(perm::Perm) = lcm(shape(perm))

"""
    sign(perm::Perm)

Returns the sign (+1 or -1) of the permutation.
"""
sign(perm::Perm) = (-1)^(degree(perm) - length(cycles(perm)))

##  largest moved point, is identity?
"""
    last_moved(perm::Perm)

Returns the largest index moved by the permutation (0 if identity).
"""
function last_moved(perm::Perm)
    for i in reverse(eachindex(perm.list))
        perm.list[i] == i || return i
    end
    return 0
end

"""
    isidentity(perm::Perm)

Check whether the permutation is the identity.
"""
isidentity(perm::Perm) = last_moved(perm) == 0

##  Cycles ('DFS' version:  each node i has two children: i^perm and i+1)
"""
    cycles(perm::Perm)

Decompose the permutation into disjoint cycles.
"""
function cycles(perm::Perm, x = 1, open = trues(degree(perm)), cycle = Int[], ccc = Array{Int}[])
    isempty(open) && return ccc  # all points have been seen
    if open[x]
        open[x] = false
        push!(cycle, x) # add x to current cycle
        cycles(perm, x^perm, open, cycle, ccc)
    else  # cycle finished, add to ccc if not empty
        isempty(cycle) || push!(ccc, cycle)
        x < length(open) && cycles(perm, x+1, open, [], ccc)
    end
    return ccc
end

##  transpositions
"""
    transpositions(n)

Return the list of adjacent transpositions `(j-1 j)` for j = 2 to n.
"""
transpositions(n::Int) = [Perm(n, [[j-1, j]]) for j in 2:n]

##  Shuffle:  Fisher-Yates
"""
    shuffle!(perm::Perm)

Apply the Fisher-Yates shuffle to `perm` in-place.
"""
function shuffle!(perm::Perm)
    n = degree(perm)
    for k in 1:n-1
        l = rand(k:n);
        perm.list[[l, k]] = perm.list[[k, l]]
    end
    return perm
end

"""
    rand(Perm, n)

Generate a random permutation of size `n`.
"""
rand(T::Type{X}, n::Int) where X <: Perm = shuffle!(one(Perm, n))

## Polya action
"""
    permuted(list, perm)

Apply the inverse of `perm` to a list — Polya-style.
"""
permuted(list::Array, perm::Perm) = list[inv(perm).list]

end # module
