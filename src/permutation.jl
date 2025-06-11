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
export permuted, transposition, transpositions, shuffle!

##  Perm data type
struct Perm
    list::Vector{Int}
end

## Identity constructors
Perm(n::Int) = Perm(collect(1:n))
one(Perm, n::Int) = Perm(n)
one(perm::Perm) = Perm(domain(perm))

## Basic attributes
degree(perm::Perm) = length(perm.list)
domain(perm::Perm) = 1:degree(perm)

## Equality, hash, comparison
==(perm::Perm, other::Perm) = perm.list == other.list
hash(perm::Perm, h::UInt) = hash(perm.list, h)
isless(perm::Perm, other::Perm) = perm.list < other.list

## Action (on points)
^(x::Int, perm::Perm) = perm.list[x]
^(xs::Vector{Int}, perm::Perm) = perm.list[xs]

## inverse permutation
function inv(perm::Perm)
    other = similar(perm.list)
    other[perm.list] = eachindex(other)
    return Perm(other)
end

## Group operations
*(perm::Perm, other::Perm) = Perm(perm.list^other)
/(perm::Perm, other::Perm) = perm * inv(other)
^(perm::Perm, other::Perm) = inv(other) * perm * other
function ^(perm::Perm, n::Int)
    n == 1 && return perm
    n == 0 && return one(perm)
    n < 0 && return inv(perm)^(-n)
    q, r = divrem(n, 2)
    return perm^q * perm^q * perm^r
end

## shape, order, sign
shape(perm::Perm) = sort(length.(cycles(perm)), rev=true)
order(perm::Perm) = lcm(shape(perm))
sign(perm::Perm) = (-1)^(degree(perm) - length(cycles(perm)))

## largest moved point, is identity?
function last_moved(perm::Perm)
    for i in reverse(eachindex(perm.list))
        perm.list[i] == i || return i
    end
    return 0
end
isidentity(perm::Perm) = last_moved(perm) == 0

## Cycles ('DFS' version:  each node i has two children: i^perm and i+1)
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

function Perm(n::Int, cycles::Array)
    perm = Perm(n)
    for c in cycles
        perm.list[c] = c[(1:end) .% end .+ 1]
    end
    return perm
end

## transpositions
function transposition(n::Int, j::Int, k::Int)
    tran = one(Perm, n)
    tran.list[[j, k]] = tran.list[[k, j]]
    return tran
end
transpositions(n::Int) = [transposition(n, j-1, j) for j in 2:n]

## shuffle:  Fisher-Yates
function shuffle!(perm::Perm)
    n = degree(perm)
    for k in 1:n-1
        l = rand(k:n);
        perm.list[[l, k]] = perm.list[[k, l]]
    end
    return perm
end
rand(T::Type{X}, n::Int) where X <: Perm = shuffle!(one(Perm, n))

## Polya action
permuted(list::Array, perm::Perm) = list[inv(perm).list]

end # module
