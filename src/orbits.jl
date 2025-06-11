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
export orbit_with_edges, orbit_with_images

## orbit
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

## on points
onPoints(x, a) = x^a

## on right
onRight(x, a) = x * a

## on pairs
onPairs(pair::Pair, a::Perm) = Pair(pair.first^a, pair.second^a)

## on sets
onSets(set::Set, a ::Perm) = Set([x^a for x in set])

##  Orbit data type
struct Orbit
    group
    elts # assume sorted!
end

## size
size(o::Orbit) = length(o.elts)

## membership, equality, comparison
in(x, o::Orbit) = x in o.elts
==(o::Orbit, other::Orbit) = o.elts[1] == other.elts[1]
isless(o::Orbit, other::Orbit) = o.elts[1] < other.elts[1]

## orbit with words
onWords(word, s) = vcat(word, s)

function orbit_with_words(aaa, x, under)
    list = [x]
    words = Array{Int}[[]]
    for (i, y) in enumerate(list)
        for (k, a) in enumerate(aaa)
            z = under(y, a)
            w = onWords(words[i], k)
            z in list || (push!(list, z), push!(words, w))
        end
    end
    return (list = list, words = words)
end

## orbitx
function orbitx(aaa, xxx, under)
    list = copy(xxx)
    for y in list
        for a in aaa
            z = under(y, a)
            z in list || push!(list, z)
        end
    end
    return list
end

## orbit with transversal
function orbit_with_transversal(aaa, x, under)
    list = [x]
    reps = [aaa[1]^0]
    for (i, y) in enumerate(list)
        for a in aaa
            z = under(y, a)
            t = reps[i] * a
            z in list || (push!(list, z), push!(reps, t))
        end
    end
    return (list = list, reps = reps)
end

## orbit with stabilizer
function orbit_with_stabilizer(aaa, x, under)
    list = [x]
    reps = [aaa[1]^0]
    stab = []
    for (i, y) in enumerate(list)
        for a in aaa
            z = under(y, a)
            t = reps[i] * a
            l = findfirst(==(z), list)
            if isnothing(l)
                push!(list, z)
                push!(reps, t)
            else   # x^(reps[i] * a) = x^reps[l]
                push!(stab, t / reps[l])
            end
        end
    end
    return (list = list, reps = reps, stab = stab)
end

## orbit with edges
function orbit_with_edges(aaa, x, under)
    list = [x]
    edges = []
    for (i, y) in enumerate(list)
        for a in aaa
            z = under(y, a)
            l = findfirst(==(z), list)
            if isnothing(l)
                push!(list, z)
                l = length(list)
            end
            push!(edges, (i, l))
        end
    end
    return (list = list, edges = edges)
end

## orbit with images
function orbit_with_images(aaa, x, under)
    list = [x]
    images = [Int[] for a in aaa]
    for (i, y) in enumerate(list)
        for (k, a) in enumerate(aaa)
            z = under(y, a)
            l = findfirst(==(z), list)
            if l == nothing
                push!(list, z);
                l = length(list);
            end
            push!(images[k], l)
        end
    end
    return (list = list, images = images)
end

end # module
