#############################################################################
##
#A  orbits.jl                                                         OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  Orbit algorithms, standard actions
##
module orbits

using ..permutation

import Base: in, isless, size, ==

export Orbit
export orbit, onPoints, onRight, onWords, onPairs, onSets
export orbit_with_words, orbit_with_transversal, orbit_with_stabilizer
export orbit_with_dist, orbit_with_tree, orbit_with_edges, orbit_with_images
export orbitl, orbitx, orbitx_with_words, orbitx_with_dist, orbitx_with_edges
export edges_from_images

## orbit
"""
    orbit(aaa, x, under)

Compute the orbit of a point `x` under the set of operators `aaa`
using the action function `under(y, a)`, which defines how elements of
the domain are acted upon.

Returns the list of points reachable from `x` by repeatedly applying
the elements of `aaa`.

# Examples
```jldoctest
julia> using OrbitAl

julia> s = Perm([2,1,3]); t = Perm([1,3,2]);

julia> orbit([s, t], 1, onPoints)
3-element Vector{Int64}:
 1
 2
 3
```
"""
function orbit(aaa, x, under::Function)
    list = [x]
    for y in list
        for a in aaa
            z = under(y, a)
            z in list || push!(list, z)
        end
    end
    return list
end

"""
    orbitl(aaa, x, under)

Like `orbit`, but uses a `Dict` for O(1) membership testing instead of
linear search. Prints a dot to stdout for each new element discovered.
Returns the list of orbit elements.

# Examples
```jldoctest
julia> using OrbitAl

julia> orbitl([Perm([2,1])], 1, onPoints)
.2-element Vector{Int64}:
 1
 2
```
"""
function orbitl(aaa, x, under::Function)
    list = [x]
    index = Dict(x => 1)
    for y in list
        for a in aaa
            z = under(y, a)
            get!(index, z) do
                push!(list, z)
                print(".")
                length(list)
            end
        end
    end
    return list
end


## common actions (to be used as action function `under`):

"""
    onPoints(x, a)

Act on a point `x` by a permutation `a`, returning `x^a`.

# Examples
```jldoctest
julia> using OrbitAl

julia> onPoints(2, Perm([2,1,3]))
1
```
"""
onPoints(x, a) = x^a

"""
    onRight(x, a)

Act on an element `x` by right-multiplying by `a`, returning `x * a`.

# Examples
```jldoctest
julia> using OrbitAl

julia> s = Perm([2,1,3]); t = Perm([1,3,2]);

julia> onRight(s, t)
Perm([3, 1, 2])
```
"""
onRight(x, a) = x * a

"""
    onPairs(pair, a)

Act on a `Pair` by applying `a` to both components: `first^a => second^a`.

# Examples
```jldoctest
julia> using OrbitAl

julia> onPairs(1 => 2, Perm([2,1,3]))
2 => 1
```
"""
onPairs(pair::Pair, a) = pair.first^a => pair.second^a

"""
    onSets(set, a)

Act on a `Set` by applying `a` to each element: `{x^a : x in set}`.

# Examples
```jldoctest
julia> using OrbitAl

julia> sort(collect(onSets(Set([1, 3]), Perm([2, 1, 3]))))
2-element Vector{Int64}:
 2
 3
```
"""
onSets(set::Set, a) = Set([x^a for x in set])

"""
    onWords(word, s)

Extend a word (index sequence) by appending generator index `s`.

# Examples
```jldoctest
julia> using OrbitAl

julia> onWords([1, 2], 1)
3-element Vector{Int64}:
 1
 2
 1
```
"""
onWords(word::Vector, s) = [word; s]

## orbit with words
"""
    orbit_with_words(aaa, x, under)

Compute the orbit of `x` under the action of the generators `aaa`,
using the action function `under` and recording the sequence
of generator indices that led to each orbit element.

Returns a named tuple:
- `list`: orbit elements
- `words`: list of index sequences tracing how each was reached

# Examples
```jldoctest
julia> using OrbitAl

julia> s = Perm([2,1,3]); t = Perm([1,3,2]);

julia> o = orbit_with_words([s, t], 1, onPoints);

julia> o.list
3-element Vector{Int64}:
 1
 2
 3

julia> o.words
3-element Vector{Vector{Int64}}:
 []
 [1]
 [1, 2]
```
"""
function orbit_with_words(aaa, x, under::Function)
    list = [x]
    words = [Int[]] # empty word Int[] maps x to x
    for (i, y) in enumerate(list)
        for (k, a) in enumerate(aaa)
            z = under(y, a)
            w = onWords(words[i], k)
            z in list || begin
                push!(list, z)
                push!(words, w)
            end
        end
    end
    return (list = list, words = words)
end

"""
    orbit_with_dist(aaa, x, under)

Compute the orbit of `x` under the generators `aaa` and record the
BFS distance from `x` to each orbit element. Returns a named tuple:
- `list`: the orbit elements
- `dist`: the corresponding distances (word lengths) from `x`

# Examples
```jldoctest
julia> using OrbitAl

julia> s = Perm([2,1,3]); t = Perm([1,3,2]);

julia> o = orbit_with_dist([s, t], 1, onPoints);

julia> o.dist
3-element Vector{Int64}:
 0
 1
 2
```
"""
function orbit_with_dist(aaa, x, under::Function)
    list = [x]
    dist = [0]
    for (i, y) in enumerate(list)
        for a in aaa
            z = under(y, a)
            z in list || begin
                push!(list, z)
                push!(dist, dist[i]+1)
            end
        end
    end
    return (list = list, dist = dist)
end

## orbit with Schreier tree
"""
    orbit_with_tree(aaa, x, under)

Compute the Schreier tree of the orbit of `x` under the generators `aaa`
using the action function `under`.

Returns a named tuple:
- `list`: a list of orbit elements
- `tree`: a `Dict` mapping each non-root node `z` to `(k => y)`, where:
    - `k` is the generator index used to reach `z`
    - `y` is the parent node so that `z = under(y, aaa[k])`

# Examples
```jldoctest
julia> using OrbitAl

julia> s = Perm([2,1,3]); t = Perm([1,3,2]);

julia> o = orbit_with_tree([s, t], 1, onPoints);

julia> o.tree[2]
1 => 1

julia> o.tree[3]
2 => 2
```
"""
function orbit_with_tree(aaa, x, under::Function)
    list = [x]
    tree = Dict()
    for (i, y) in enumerate(list)
        for (k, a) in enumerate(aaa)
            z = under(y, a)
            z in list || begin
                push!(list, z)
                tree[z] = k => y  # combined label and parent
            end
        end
    end
    return (list = list, tree = tree)
end

## orbit with transversal
"""
    orbit_with_transversal(aaa, x, under)

Compute the orbit of `x` under the action of `aaa`, along with
coset representatives. Returns a named tuple:
- `list`: the orbit of `x`
- `reps`: permutations such that `x^reps[i] == list[i]` for each `i`

# Examples
```jldoctest
julia> using OrbitAl

julia> s = Perm([2,1,3]); t = Perm([1,3,2]);

julia> o = orbit_with_transversal([s, t], 1, onPoints);

julia> [1^r for r in o.reps]
3-element Vector{Int64}:
 1
 2
 3
```
"""
function orbit_with_transversal(aaa, x, under::Function)
    list = [x]
    reps = [aaa[1]^0] # identity maps x to x
    for (i, y) in enumerate(list)
        for a in aaa
            z = under(y, a)
            t = onRight(reps[i], a)
            z in list || begin
                push!(list, z)
                push!(reps, t)
            end
        end
    end
    return (list = list, reps = reps)
end

## orbit with stabilizer
"""
    orbit_with_stabilizer(aaa, x, under)

Compute the orbit of `x` under generators `aaa` using Schreier's Lemma.
Returns a named tuple:
- `list`: the orbit elements
- `reps`: permutations such that `x^reps[i] == list[i]`
- `stab`: Schreier generators for the stabilizer of `x` (may contain duplicates)

# Examples
```jldoctest
julia> using OrbitAl

julia> s = Perm([2,1,3]); t = Perm([1,3,2]);

julia> o = orbit_with_stabilizer([s, t], 1, onPoints);

julia> length(o.list)
3

julia> unique(o.stab)
1-element Vector{Perm}:
 Perm([1, 3, 2])
```
"""
function orbit_with_stabilizer(aaa, x, under::Function)
    list = [x]
    index = Dict(x => 1)
    reps = [aaa[1]^0]  # identity maps x to x
    stab = Perm[]
    for (i, y) in enumerate(list)
        for a in aaa
            z = under(y, a)
            t = onRight(reps[i], a)
            l = get!(index, z) do
                push!(list, z)
                push!(reps, t)
                length(list)
            end   # x^(reps[i] * a) = x^reps[l]
            t == reps[l] || push!(stab, t / reps[l])
        end
    end
    return (list = list, reps = reps, stab = stab)
end

## orbit with edges
"""
    orbit_with_edges(aaa, x, under)

Compute the orbit of `x` under the generators `aaa` and return:
- `list`: the list of orbit elements
- `edges`: directed edge pairs `(i, j)` (both directions included for each
  generator application `list[i] -> list[j]` with `i ≠ j`)

# Examples
```jldoctest
julia> using OrbitAl

julia> o = orbit_with_edges([Perm([2,1])], 1, onPoints);

julia> o.list
2-element Vector{Int64}:
 1
 2

julia> sort(o.edges)
2-element Vector{Any}:
 (1, 2)
 (2, 1)
```
"""
function orbit_with_edges(aaa, x, under::Function)
    list = [x]
    index = Dict(x => 1)
    edges = Set()
    for (i, y) in enumerate(list)
        for a in aaa
            z = under(y, a)
            l = get!(index, z) do
                push!(list, z)
                length(list)
            end
            i == l || push!(edges, (i, l))
        end
    end
    return (list = list, edges = collect(edges))
end

"""
    orbit_with_images(aaa, x, under)

Compute the orbit of `x` under the generators `aaa` and record, for
each generator, the list of image indices. Returns a named tuple:
- `list`: the orbit elements
- `images`: a vector of index lists, one per generator, where `images[k][i]`
  is the index of `under(list[i], aaa[k])` in `list`

# Examples
```jldoctest
julia> using OrbitAl

julia> s = Perm([2,1,3]); t = Perm([1,3,2]);

julia> o = orbit_with_images([s, t], 1, onPoints);

julia> o.images[1]
3-element Vector{Int64}:
 2
 1
 3

julia> o.images[2]
3-element Vector{Int64}:
 1
 3
 2
```
"""
function orbit_with_images(aaa, x, under::Function)
    list = [x]
    index = Dict(x => 1)
    images = [Int[] for a in aaa]
    for (i, y) in enumerate(list)
        for (k, a) in enumerate(aaa)
            z = under(y, a)
            l = get!(index, z) do
                push!(list, z)
                length(list)
            end
            push!(images[k], l)
        end
    end
    return (list = list, images = images)
end

## edges from images
"""
    edges_from_images(images)

Convert the `images` field from [`orbit_with_images`](@ref) into a list of
labeled directed edges. Returns a vector of triples `(i, j, k)` where:
- `i` is the source index
- `j` is the target index (`images[k][i] == j`)
- `k` is the generator index

Self-loops (`i == j`) are excluded.

# Examples
```jldoctest
julia> using OrbitAl

julia> s = Perm([2,1,3]); t = Perm([1,3,2]);

julia> o = orbit_with_images([s, t], 1, onPoints);

julia> edges_from_images(o.images)
4-element Vector{Tuple{Int64, Int64, Int64}}:
 (1, 2, 1)
 (2, 1, 1)
 (2, 3, 2)
 (3, 2, 2)
```
"""
function edges_from_images(images)
    return [
        (i, j, k) for (k,img) in enumerate(images)
        for (i,j) in enumerate(img) if i != j
    ]
end

"""
    orbitx(aaa, xxx, under)

Like `orbit`, but starts from a set of points `xxx` rather than a single
point. Returns the list of all elements reachable from any element of `xxx`.

# Examples
```jldoctest
julia> using OrbitAl

julia> s = Perm([2,1,3]); t = Perm([1,3,2]);

julia> orbitx([s, t], [1, 2], onPoints)
3-element Vector{Int64}:
 1
 2
 3
```
"""
function orbitx(aaa, xxx, under::Function)
    list = xxx
    for y in list
        for a in aaa
            z = under(y, a)
            z in list || push!(list, z)
        end
    end
    return list
end

"""
    orbitx_with_words(aaa, xxx, under)

Like `orbit_with_words`, but starts from a set of points `xxx`. Each
starting point `xxx[i]` gets initial word `[i]`, and subsequent
generator applications are appended. Returns a named tuple:
- `list`: all reachable elements
- `words`: index sequences where `words[j][1]` is the seed index into `xxx`
  and `words[j][2:end]` are the generator indices applied to reach `list[j]`

# Examples
```jldoctest
julia> using OrbitAl

julia> s = Perm([2,1,3]); t = Perm([1,3,2]);

julia> o = orbitx_with_words([s, t], [1, 2], onPoints);

julia> o.list
3-element Vector{Int64}:
 1
 2
 3

julia> o.words
3-element Vector{Vector{Int64}}:
 [1]
 [2]
 [2, 2]
```
"""
function orbitx_with_words(aaa, xxx, under::Function)
    list = xxx
    words = [[i] for i in eachindex(xxx)]
    index = Dict(x => i for (i, x) in enumerate(xxx))
    for (i, y) in enumerate(list)
        for (k, a) in enumerate(aaa)
            z = under(y, a)
            w = onWords(words[i], k)
            get!(index, z) do
                push!(list, z)
                push!(words, w)
                length(list)
            end
        end
    end
    return (list = list, words = words)
end

function orbitx_with_dist(aaa, xxx, under::Function)
    list = xxx
    dist = [0 for x in xxx]
    for (i, y) in enumerate(list)
        for a in aaa
            z = under(y, a)
            z in list || begin
                push!(list, z)
                push!(dist, dist[i]+1)
            end
        end
    end
    return (list = list, dist = dist)
end



"""
    orbitx_with_edges(aaa, xxx, under)

Like `orbit_with_edges`, but starts from a set of points `xxx`. Returns
a named tuple:
- `list`: all reachable elements
- `edges`: directed edge pairs `(i, j)` (both directions included)

# Examples
```jldoctest
julia> using OrbitAl

julia> s = Perm([2,1,3]); t = Perm([1,3,2]);

julia> o = orbitx_with_edges([s, t], [1, 2], onPoints);

julia> o.list
3-element Vector{Int64}:
 1
 2
 3

julia> sort(o.edges)
4-element Vector{Any}:
 (1, 2)
 (2, 1)
 (2, 3)
 (3, 2)
```
"""
function orbitx_with_edges(aaa, xxx, under::Function)
    list = xxx
    index = Dict(x => i for (i, x) in enumerate(xxx))
    edges = Set()
    for (i, y) in enumerate(list)
        for (k, a) in enumerate(aaa)
            z = under(y, a)
            l = get!(index, z) do
                push!(list, z)
                length(list)
            end
            i == l || push!(edges, (i, l))
        end
    end
    return (list = list, edges = collect(edges))
end


"""
    Orbit

Represents the conjugacy class (or similar orbit) of an element under a
group action. Created by `perm ^ group`, which computes the orbit of `perm`
under conjugation. Fields:
- `group`: the acting group
- `elts`: a sorted list of orbit elements

Supports `in`, `==` (by minimum element), `isless` (by minimum element),
and `size`.

# Examples
```jldoctest
julia> using OrbitAl

julia> s = Perm([2,1,3]); t = Perm([1,3,2]);

julia> G = PermGp([s, t], one(s));

julia> cc = conjClasses(G);

julia> length(cc)
3

julia> cc[1].elts
1-element Vector{Perm}:
 Perm([1, 2, 3])
```
"""
struct Orbit
    group
    elts # assume sorted!
end

## size, membership, equality, comparison
size(o::Orbit) = length(o.elts)
in(x, o::Orbit) = x in o.elts
==(o::Orbit, other::Orbit) = o.elts[1] == other.elts[1]
isless(o::Orbit, other::Orbit) = o.elts[1] < other.elts[1]

end # module
