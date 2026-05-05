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
export orbitl, orbitx, orbitx_with_words, orbitx_with_edges

## orbit
"""
    orbit(aaa, x, under)

Compute the orbit of a point `x` under the set of operators `aaa`
using the action function `under(y, a)`, which defines how elements of
the domain are acted upon.

Returns the list of points reachable from `x` by repeatedly applying
the elements of `aaa`.
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

Like `orbit`, but prints a dot for each new element discovered.
Returns the list of orbit elements.
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
"""
onPoints(x, a) = x^a

"""
    onRight(x, a)

Act on an element `x` by right-multiplying by `a`, returning `x * a`.
"""
onRight(x, a) = x * a

"""
    onPairs(pair, a)

Act on a `Pair` by applying `a` to both components: `first^a => second^a`.
"""
onPairs(pair::Pair, a) = pair.first^a => pair.second^a

"""
    onSets(set, a)

Act on a `Set` by applying `a` to each element: `{x^a : x in set}`.
"""
onSets(set::Set, a) = Set([x^a for x in set])

"""
    onWords(word, s)

Extend a word (index sequence) by appending generator index `s`.
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

Returns:
- `list`: a list of orbit elements
- `tree`: a `Dict` mapping each non-root node `z` to `(k => y)`, where:
    - `k` is the generator index used to reach `z`
    - `y` is the parent node so that `z = under(y, aaa[k])`
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
coset representatives (as permutations). Returns a named tuple
with fields:
- `list`: the orbit of `x`
- `reps`: corresponding permutations mapping `x` to each orbit element.
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

Compute the orbit of `x` under generators `aaa`, recording:
- `list`: the orbit elements,
- `reps`: permutations mapping `x` to each orbit element,
- `stab`: generators of the stabilizer of `x` (via Schreier's Lemma).

Returns a named tuple `(list, reps, stab)`.
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
- `edges`: a list of undirected edges (as sorted index pairs)
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
starting point `xxx[i]` is assigned the initial word `[i]`. Returns a
named tuple:
- `list`: all reachable elements
- `words`: the index sequence tracing how each was reached
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

"""
    orbitx_with_edges(aaa, xxx, under)

Like `orbit_with_edges`, but starts from a set of points `xxx`. Returns
a named tuple:
- `list`: all reachable elements
- `edges`: undirected edges as sorted index pairs
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

Represents the orbit of an element under a group action. Fields:
- `group`: the group acting
- `elts`: a sorted list of orbit elements

Supports `in`, `==`, `isless` (by first element), and `size`.
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
