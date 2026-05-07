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

## [l_1, l_2 ..  l_n] -> [l_1, l_2 - l_1, l_3 - l_2 .. l_n - l_{n-1}]
"""
    newtonDif(list)

Compute the Newton forward differences of `list`. Returns a vector whose first
element is `list[1]` and whose remaining elements are the differences between
consecutive elements.

Inverse of `newtonSum`: `newtonDif(newtonSum(list)) == list`.

# Examples
```jldoctest
julia> using OrbitAl

julia> newtonDif([1, 3, 6, 10])
4-element Vector{Int64}:
 1
 2
 3
 4
```
"""
function newtonDif(list::Vector)
    isempty(list) && return []
    return [list[1]; diff(list)]
end

## [s_1, s_2 ..  s_n] -> [s_1, s_1 + s_2, s_1 + s_2 + s_3 .. s_1 + ... + s_n]
"""
    newtonSum(list)

Compute the Newton partial (cumulative) sums of `list`. The `k`-th element of
the result is `list[1] + list[2] + … + list[k]`.

Inverse of `newtonDif`: `newtonSum(newtonDif(list)) == list`.

# Examples
```jldoctest
julia> using OrbitAl

julia> newtonSum([1, 2, 3, 4])
4-element Vector{Int64}:
  1
  3
  6
 10
```
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

Right-handed Newton differences: the mirror image of `newtonDif` applied
right-to-left. For a non-increasing sequence (such as a partition), the `k`-th
element of the result is `list[k] - list[k+1]` (with `list[end]` as the last
element), giving a non-negative vector.

Inverse of `newtonSumR`: `newtonDifR(newtonSumR(list)) == list`.

# Examples
```jldoctest
julia> using OrbitAl

julia> newtonDifR([6, 5, 3])
3-element Vector{Int64}:
 1
 2
 3
```
"""
newtonDifR(list::Vector) = reverse(newtonDif(reverse(list)))

## [s_1, s_2 ..  s_n] -> [s_1 + ... + s_n, s_2 + ... + s_n .. s_n]
"""
    newtonSumR(list)

Right-handed Newton partial sums: the mirror image of `newtonSum` applied
right-to-left. The `k`-th element of the result is
`list[k] + list[k+1] + … + list[end]`.

Inverse of `newtonDifR`: `newtonSumR(newtonDifR(list)) == list`.

# Examples
```jldoctest
julia> using OrbitAl

julia> newtonSumR([1, 2, 3])
3-element Vector{Int64}:
 6
 5
 3
```
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
    compositionSubset(n, set)

Given a subset `set ⊆ 1:n-1`, return the composition of `n` whose partial sums
form the complement of `set` in `1:n`.

Inverse of `subsetComposition`.

# Examples
```jldoctest
julia> using OrbitAl

julia> compositionSubset(9, [4, 5, 7])
6-element Vector{Int64}:
 1
 1
 1
 3
 2
 1
```
"""
compositionSubset(n::Int, set::Vector{Int}) = newtonDif(setdiff(1:n, set))

##  how to turn a composition of n into a subset of [1..n-1]
"""
    subsetComposition(n, com)

Given a composition `com` of `n`, return the subset of `1:n-1` corresponding
to the complement of the partial sums of `com` in `1:n`.

Inverse of `compositionSubset`.

# Examples
```jldoctest
julia> using OrbitAl

julia> subsetComposition(9, [1, 1, 1, 3, 2, 1])
3-element Vector{Int64}:
 4
 5
 7
```
"""
subsetComposition(n::Int, com::Vector{Int}) = setdiff(1:n, newtonSum(com))

##  examples of orbits:

##  take-away action
"""
    takeAway(set, s)

Remove element `s` from `set` and return the resulting set.

# Examples
```jldoctest
julia> using OrbitAl

julia> takeAway([1, 2, 3], 2)
2-element Vector{Int64}:
 1
 3
```
"""
takeAway(set::Union{Set,Vector}, s) = setdiff(set, s)

##  power set as orbit under take-away
"""
    subsets(set)

Compute all subsets of `set` as the orbit of `set` under the take-away action.
Returns a vector of all subsets (as vectors), starting with `set` itself.

# Examples
```jldoctest
julia> using OrbitAl

julia> subsets([1, 2])
4-element Vector{Vector{Int64}}:
 [1, 2]
 [2]
 [1]
 []
```
"""
subsets(set) = orbit(set, set, takeAway)

##  partitions of n as takeAway orbit of composition classes
"""
    partitions(n)

Enumerate all integer partitions of `n` as an orbit under a partition-removal
action. Returns a vector of non-increasing integer vectors, in reverse
lexicographic order (smallest partition first).

# Examples
```jldoctest
julia> using OrbitAl

julia> partitions(4)
5-element Vector{Vector{Int64}}:
 [1, 1, 1, 1]
 [2, 1, 1]
 [2, 2]
 [3, 1]
 [4]
```
"""
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
        if isnothing(pos)
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

"""
    standardYTs(lambda)

Return all standard Young tableaux of shape `lambda` (a partition given as a
non-increasing integer vector). Each tableau is represented as a path in the
Young lattice: a vector of `(row, col)` pairs recording where each number
`1, 2, …, n` was placed.

# Examples
```jldoctest
julia> using OrbitAl

julia> length(standardYTs([2, 1]))
2

julia> length(standardYTs([3, 2, 1]))
16
```
"""
function standardYTs(lambda)
    next = youngLattice(lambda).next
    return pathfinder(1, next, Any[missing for x in next])
end

"""
    tableau_path(path)

Convert a standard Young tableau path (a vector of `(row, col)` pairs, as
returned by `standardYTs`) into a row array, where `tab[r][c]` is the entry
placed in row `r`, column `c`.

# Examples
```jldoctest
julia> using OrbitAl

julia> paths = standardYTs([2, 1]);

julia> tableau_path(paths[1])
2-element Vector{Array{Int64}}:
 [1, 2]
 [3]

julia> tableau_path(paths[2])
2-element Vector{Array{Int64}}:
 [1, 3]
 [2]
```
"""
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
