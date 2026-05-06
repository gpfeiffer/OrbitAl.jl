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

# Examples
```jldoctest
julia> using OrbitAl

julia> p = Perm([6, 4, 7, 2, 5, 9, 8, 3, 1])
Perm([6, 4, 7, 2, 5, 9, 8, 3, 1])

julia> p == Perm(9, [[1, 6, 9], [2, 4], [3, 7, 8]])
true

julia> id9 = Perm(9)
Perm([1, 2, 3, 4, 5, 6, 7, 8, 9])

julia> p == id9
false
```
"""
struct Perm
    list::Vector{Int}
end

##  Constructors
"""
    Perm(n::Int)

Construct the identity permutation on `1:n`.

# Examples
```jldoctest
julia> using OrbitAl

julia> Perm(3)
Perm([1, 2, 3])
```
"""
Perm(n::Int) = Perm(collect(1:n))

"""
    one(Perm, n::Int)

Alias for `Perm(n)` — the identity permutation on `n` elements.

# Examples
```jldoctest
julia> using OrbitAl

julia> one(Perm, 3)
Perm([1, 2, 3])
```
"""
one(Perm, n::Int) = Perm(n)

"""
    one(perm::Perm)

Identity permutation with the same degree as `perm`.

# Examples
```jldoctest
julia> using OrbitAl

julia> one(Perm([2, 1, 3]))
Perm([1, 2, 3])
```
"""
one(perm::Perm) = Perm(domain(perm))

"""
    Perm(n::Int, cycles::Array)

Construct a permutation of degree `n` from disjoint cycles.
Each cycle is a vector of indices.

# Examples
```jldoctest
julia> using OrbitAl

julia> Perm(4, [[1, 2, 3]])
Perm([2, 3, 1, 4])
```
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

# Examples
```jldoctest
julia> using OrbitAl

julia> degree(Perm([2, 1, 3]))
3
```
"""
degree(perm::Perm) = length(perm.list)

"""
    domain(perm::Perm)

Returns the domain of the permutation as range `1:n`.

# Examples
```jldoctest
julia> using OrbitAl

julia> domain(Perm([2, 1, 3]))
1:3
```
"""
domain(perm::Perm) = 1:degree(perm)

##  Equality, hash, comparison
"""
    ==(p::Perm, q::Perm)

Check two permutations for equality.

# Examples
```jldoctest
julia> using OrbitAl

julia> Perm([2, 1, 3]) == Perm([2, 1, 3])
true

julia> Perm([2, 1, 3]) == Perm([1, 2, 3])
false
```
"""
==(perm::Perm, other::Perm) = perm.list == other.list

"""
    hash(perm::Perm, h::UInt)

Compute the hash of a permutation, consistent with `==`.
Enables use of `Perm` as a dictionary key or set element.

# Examples
```jldoctest
julia> using OrbitAl

julia> hash(Perm([2, 1, 3])) == hash(Perm([2, 1, 3]))
true
```
"""
hash(perm::Perm, h::UInt) = hash(perm.list, h)

"""
    isless(p::Perm, q::Perm)

Lexicographic comparison of two permutations (by image list).

# Examples
```jldoctest
julia> using OrbitAl

julia> isless(Perm([1, 2, 3]), Perm([2, 1, 3]))
true
```
"""
isless(perm::Perm, other::Perm) = perm.list < other.list

##  Action (on points)

"""
    x ^ p

Apply the permutation `p` to the integer point `x` or integer vector `xs`,
returning the image under `p`.

# Examples
```jldoctest
julia> using OrbitAl

julia> 1 ^ Perm([2, 1, 3])
2

julia> [1, 2, 3] ^ Perm([2, 1, 3])
3-element Vector{Int64}:
 2
 1
 3
```
"""
^(x::Int, perm::Perm) = perm.list[x]
^(xs::Vector{Int}, perm::Perm) = perm.list[xs]

##  Inverse permutation
"""
    inv(perm::Perm)

Returns the inverse of the permutation.

# Examples
```jldoctest
julia> using OrbitAl

julia> p = Perm([2, 3, 1]);

julia> inv(p)
Perm([3, 1, 2])

julia> p * inv(p) == one(p)
true
```
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

# Examples
```jldoctest
julia> using OrbitAl

julia> Perm([2, 1, 3]) * Perm([1, 3, 2])
Perm([3, 1, 2])
```
"""
*(perm::Perm, other::Perm) = Perm(perm.list^other)

"""
    p / q

Right division: `p * inv(q)`.

# Examples
```jldoctest
julia> using OrbitAl

julia> p = Perm([2, 3, 1]); q = Perm([1, 3, 2]);

julia> p / q == p * inv(q)
true
```
"""
/(perm::Perm, other::Perm) = perm * inv(other)

"""
    p ^ q

Conjugation of `p` by `q`: `inv(q) * p * q`.

# Examples
```jldoctest
julia> using OrbitAl

julia> p = Perm([2, 1, 3]); q = Perm([2, 3, 1]);

julia> p^q == inv(q) * p * q
true
```
"""
^(perm::Perm, other::Perm) = inv(other) * perm * other

"""
    p ^ n

Raise the permutation `p` to the integer power `n`.

# Examples
```jldoctest
julia> using OrbitAl

julia> p = Perm([2, 3, 1]);

julia> p^3 == one(p)
true
```
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

Returns the cycle type: cycle lengths sorted in decreasing order.

# Examples
```jldoctest
julia> using OrbitAl

julia> shape(Perm([2, 3, 1, 5, 4]))
2-element Vector{Int64}:
 3
 2
```
"""
shape(perm::Perm) = sort(length.(cycles(perm)), rev=true)

"""
    order(perm::Perm)

Returns the order of the permutation (the smallest positive `n` with `p^n == one(p)`).

# Examples
```jldoctest
julia> using OrbitAl

julia> order(Perm([2, 3, 1, 5, 4]))
6
```
"""
order(perm::Perm) = lcm(shape(perm))

"""
    sign(perm::Perm)

Returns the sign of the permutation: `+1` for even permutations, `-1` for odd.

# Examples
```jldoctest
julia> using OrbitAl

julia> sign(Perm([2, 1, 3]))
-1

julia> sign(Perm([2, 3, 1]))
1
```
"""
sign(perm::Perm) = (-1)^(degree(perm) - length(cycles(perm)))

##  largest moved point, is identity?
"""
    last_moved(perm::Perm)

Returns the largest index moved by the permutation, or `0` if it is the identity.

# Examples
```jldoctest
julia> using OrbitAl

julia> last_moved(Perm([2, 1, 3]))
2

julia> last_moved(Perm([1, 2, 3]))
0
```
"""
function last_moved(perm::Perm)
    for i in reverse(eachindex(perm.list))
        perm.list[i] == i || return i
    end
    return 0
end

"""
    isidentity(perm::Perm)

Return `true` if the permutation is the identity.

# Examples
```jldoctest
julia> using OrbitAl

julia> isidentity(Perm([1, 2, 3]))
true

julia> isidentity(Perm([2, 1, 3]))
false
```
"""
isidentity(perm::Perm) = last_moved(perm) == 0

##  Cycles ('DFS' version:  each node i has two children: i^perm and i+1)
"""
    cycles(perm::Perm)

Decompose the permutation into disjoint cycles, returned as a vector of
index vectors in order of their smallest element.

# Examples
```jldoctest
julia> using OrbitAl

julia> cycles(Perm([2, 3, 1, 5, 4]))
2-element Vector{Array{Int64}}:
 [1, 2, 3]
 [4, 5]
```
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

Return the list of adjacent transpositions `(j-1 j)` for `j = 2, …, n`.
These are the standard Coxeter generators of the symmetric group `Sₙ`.

# Examples
```jldoctest
julia> using OrbitAl

julia> transpositions(3)
2-element Vector{Perm}:
 Perm([2, 1, 3])
 Perm([1, 3, 2])
```
"""
transpositions(n::Int) = [Perm(n, [[j-1, j]]) for j in 2:n]

##  Shuffle:  Fisher-Yates
"""
    shuffle!(perm::Perm)

Apply the Fisher-Yates shuffle to `perm` in-place and return it.

# Examples
```jldoctest
julia> using OrbitAl

julia> p = Perm(5); shuffle!(p);

julia> sort(p.list)
5-element Vector{Int64}:
 1
 2
 3
 4
 5
```
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

Generate a uniformly random permutation of size `n`.

# Examples
```jldoctest
julia> using OrbitAl

julia> p = rand(Perm, 5);

julia> sort(p.list)
5-element Vector{Int64}:
 1
 2
 3
 4
 5
```
"""
rand(T::Type{X}, n::Int) where X <: Perm = shuffle!(one(Perm, n))

## Polya action
"""
    permuted(list, perm)

Return `list` rearranged by the inverse of `perm`: `permuted(list, perm)[i] == list[i^inv(perm)]`.
This is the standard Pólya action of a permutation on a list.

# Examples
```jldoctest
julia> using OrbitAl

julia> permuted([10, 20, 30], Perm([2, 3, 1]))
3-element Vector{Int64}:
 30
 10
 20
```
"""
permuted(list::Array, perm::Perm) = list[inv(perm).list]

end # module
