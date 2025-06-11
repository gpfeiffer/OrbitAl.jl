#############################################################################
##
#A  syt.jl                                                            OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  Partitions and Standard Young Tableaus as orbits under suitable actions
##
module syt

using ..orbits

export newtonDif, newtonSum, newtonDifR, newtonSumR
export compositionSubset, subsetComposition
export takeAway, subsets, partitions, standardYTs, tableau_path

#############################################################################
##
##  newtonSum(list) and newtonDif(list)
##
##  Newton's partial sums and forward differences of a list
##
##  properties:
##    newtonDif(newtonSum(list)) = list
##    newtonSum(newtonDif(list)) = list
##

## [l_1, l_2 ..  l_n] -> [l_1, l_2 - l_1, l_3 - l_2 .. l_n -l_{n-1}]
"""
    newtonDif(list)

Compute the Newton forward difference of the input list.
Returns a list of first element and the differences between subsequent
elements.
"""
function newtonDif(list::Vector)
    isempty(list) && return []
    return [list[1]; diff(list)]
end

## [s_1, s_2 ..  s_n] -> [s_1, s_1 + s_2, s_1 + s_2 + s_3 .. s_1 + ... + s_n]
"""
    newtonSum(list)

Compute the Newton partial sums of the input list.
Returns a list where each element is the cumulative sum up to that
index.
"""
function newtonSum(list::Vector)
    isempty(list) && return []
    return cumsum(list)
end

##  right handed versions:
##    newtonDifR(newtonSumR(list)) = list
##    newtonSumR(newtonDifR(list)) = list
##

## [l_1, l_2 ..  l_n] -> [l_2 - l_1 .. l_n - l_{n-1}, l_n]
"""
    newtonDifR(list)

Right-handed (backward) Newton difference:
Equivalent to reversing the result of forward difference on reversed list.
"""
newtonDifR(list::Vector) = reverse(newtonDif(reverse(list)))

## [s_1, s_2 ..  s_n] -> [s_1 + ... + s_n, s_2 + ... + s_n .. s_n]
"""
    newtonSumR(list)

Right-handed (backward) Newton sum:
Equivalent to reversing the result of forward sum on reversed list.
"""
newtonSumR(list::Vector) = reverse(newtonSum(reverse(list)))

#############################################################################
##
##  Compositions (of n) vs Subsets (of {1..n-1})
##
##  how to turn a subset of [1..n-1] into a composition of n:
##  * find the complement cmp of L in [1..n]
##  * compute Newton differences of cmp
##  e.g. n = 9, L =    45 7  \subseteq [1..8],
##            cmp = 123  6 89,
##            com = 111  3 21.
##
"""
    compositionSubset(n::Int, set::Vector{Int})

Given a subset `set ⊆ 1:n-1`, return the composition of `n` determined by
the **complement** of `set` in `1:n`, interpreted via Newton differences.

Example:
    n = 9, set = [4,5,7]
    → complement = [1,2,3,6,8,9]
    → composition = [1,1,1,3,2,1]
"""
compositionSubset(n::Int, set::Vector{Int}) = newtonDif(setdiff(1:n, set))

##  how to turn a composition of n into a subset of [1..n-1]
"""
    subsetComposition(n::Int, com::Vector{Int})

Given a composition `com` of `n`, return the subset of `1:n-1` corresponding
to the **complement** of its partial sums.

This undoes `compositionSubset`.
"""
subsetComposition(n::Int, com::Vector{Int}) = setdiff(1:n, newtonSum(com))

##  examples of orbits:

##  take-away action
"""
    takeAway(set, s)

The take-away action: remove `s` from `set`, returning a new set.
"""
takeAway(set::Union{Set,Vector}, s) = setdiff(set, s)

##  power set as orbit under take-away
subsets(set) = orbit(set, set, takeAway)

##  partitions of n as takeAway orbit of composition classes
function partitions(n::Int)
    function takeAwayPartition(com::Vector{Int}, a::Int)
        com = compositionSubset(n, takeAway(subsetComposition(n, com), a))
        return sort(com, rev=true)
    end
    reverse(orbit(1:n-1, [n], takeAwayPartition))  # canonical order
end

#############################################################################
##
##  SYT: a standard Young tableau is a shortest path in the Young lattice
##

##  normal case action
function remove1Hook2(lambda, k)
    new = copy(lambda)
    new[k] -= 1
    return new
end

##  special case action
function remove1Hook1(lambda, k)
    lambda[k] == 1 && return lambda[1:k-1]
    remove1Hook2(lambda, k)
end

##  Young lattice and shortest paths.  2-step process.
##  1. build graph "bottom up" from lambda; for each edge, record the row/col
##     positions it corresponds to.
##  2. find shortest paths "top down" from 0
##  (consider setting this up as starting with a list of partitions of n)
##
function youngLattice(lambda)

    ## how to add to the list
    function grow(list, y, r, c)
        pos = findfirst(==(y), list)
        if pos == nothing
            push!(list, y)
            push!(next, [])
            pos = length(list)
        end
        push!(next[i], (pos = pos, row = r, col = c))
    end

    ## orbit-with-edges
    list = [lambda]
    next = [[] for x in list]
    i = 0
    while i < length(list)
        i += 1
        x = list[i]
        l = length(x)
        if l > 0
            dif = newtonDifR(x)   # ;-)
            grow(list, remove1Hook1(x, l), l, x[l])  # dif[l] > 0
            for k in reverse(1:l-1)
                dif[k] > 0 && grow(list, remove1Hook2(x, k), k, x[k])
            end
        end
    end
    return (list = list, next = next)
end

##  DFS with a visitor that makes shortest paths if unknown.
function pathfinder(x, next, path)
    if ismissing(path[x])
        path[x] = []
        for z in next[x]
            for p in pathfinder(z.pos, next, path)
                push!(path[x], vcat(p, [(z.row, z.col)]))
            end
        end
        path[x] == [] && push!(path[x], [])   # empty path
    end
    return path[x]
end

## list of SYTs of shape lambda
function standardYTs(lambda)
    next = youngLattice(lambda).next
    return pathfinder(1, next, Any[missing for x in next])
end

##  tableau from path
function tableau_path(path)
    tab = Array{Int}[]
    for i in 1:length(path)
        r, c = path[i]
        c == 1 && push!(tab, [])
        push!(tab[r], i)
    end
    return tab
end

end # module
