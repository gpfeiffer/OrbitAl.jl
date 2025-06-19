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
export orbit, orbitx, onPoints, onRight, onWords, onPairs, onSets
export orbit_with_words, orbit_with_transversal, orbit_with_stabilizer
export orbit_with_dist, orbit_with_tree, orbit_with_edges, orbit_with_images

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

## orbitx: multiple starting points
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

## common actions (to be used as action function `under`):

# on points
onPoints(x, a::Perm) = x^a

# on right
onRight(x, a::Perm) = x * a

# on pairs
onPairs(pair::Pair, a::Perm) = pair.first^a => pair.second^a

# on sets
onSets(set::Set, a ::Perm) = Set([x^a for x in set])

# on words
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
    words = Array{Int}[[]]  # word[i] gives path to list[i]
    for (i, y) in enumerate(list)
        for (k, a) in enumerate(aaa)
            z = under(y, a)
            w = onWords(words[i], k)
            z in list || begin
                push!(list, z); push!(words, w)
            end
        end
    end
    return (list = list, words = words)
end

## orbit with distance (= word length)
function orbit_with_dist(aaa, x, under::Function)
    list = [x]
    dist = [0]
    for (i, y) in enumerate(list)
        for (k, a) in enumerate(aaa)
            z = under(y, a)
            z in list || begin
                push!(list, z); push!(dist, dist[i]+1)
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
- `orbit`: a list of orbit elements
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
    reps = [one(aaa[1])]  # identity permutation
    for (i, y) in enumerate(list)
        for a in aaa
            z = under(y, a)
            t = onRight(reps[i], a)
            z in list || begin
                push!(list, z); push!(reps, t)
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
    reps = [one(aaa[1])]  # identity permutation
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
            i == l || push!(edges, Tuple(sort([i, l])))
        end
    end
    return (list = list, edges = collect(edges))
end

## orbit with images
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
      (i, j, k) for (k,img) in enumerate(images) for (i,j) in enumerate(img)
    ]
end

##  Orbit data type
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
