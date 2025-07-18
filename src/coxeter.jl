#############################################################################
##
#A  coxeter.jl                                                        OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  Constructing finite Coxeter groups and related orbits.
##
module coxeter

using ..permutation
using ..orbits
using ..permgroup

import Base: ^, size
import ..permutation: Perm
import ..permutation: isidentity, last_moved
import ..permgroup: PermGp

export coxeterGraph, cartanMat, CoxeterGp, coxeterConjugacyClasses
export coxeterLength, coxeterWord, permCoxeterWord, reflections
export size, prefixes, prefixes_with_edges, shape, shapes
export longestElt, parabolicTransversal

function coxeterGraph(series::String, rank::Int)
    edges = [(j-1,j) for j in 2:rank]     # type A: chain
    series >= "D" && (edges[1] = (1,3))   # type D: add fork
    series >= "E" && (edges[2] = (2,4))   # type E: another fork
    return edges
end

function cartanMat(series::String, rank::Int)
    C = [i == j ? 2 : 0 for i in 1:rank, j in 1:rank]
    for (i,j) in coxeterGraph(series, rank)
        C[i,j] = -1;  C[j,i] = -1
    end
    series == "B" && (C[1,2] = -2)
    series == "C" && (C[2,1] = -2)
    series == "F" && (C[3,4] = -2)  # sic!
    return C
end

absRoot(r) = sum(r) < 0 ? -r : r
onRoots(x, a) = absRoot(onRight(x, a))

Perm(a, xxx, under) = Perm(indexin([under(x, a) for x in xxx], xxx))

abstract type ACoxeterGp <: APermGp end

struct CoxeterGp <: ACoxeterGp
    gens::Vector{Perm}
    one::Perm
    data::Dict{Symbol, Any}
end

function reflection(C::Matrix{Int}, s::Int)
    mat = C^0
    mat[:,s] = mat[s,:] - C[s,:]
    return mat
end

function CoxeterGp(C::Matrix{Int})
    one = C^0
    S = axes(C, 1)
    mats = [reflection(C, s) for s in S]
    roots = orbitx_with_words(mats, [one[i:i,:] for i in S], onRoots)
    data = Dict(:mats => mats, :roots => roots, :rank => length(S))
    data[:C] = C;  data[:N] = length(roots.list)
    data[:phi] = [roots.list; -roots.list]
    perms = [Perm(m, data[:phi], onRight) for m in mats]
    return CoxeterGp(perms, perms[1]^0, data)
end

function reflections(W)
    refl(w) = W.gens[w[1]]^prod(W.gens[w[2:end]]; init=W.one)
    refl.(W.data[:roots].words)
end

function coxeterLength(W, w)
    N = W.data[:N]
    return count(i^w > N for i in 1:N)
end

permCoxeterWord(W, word) = prod(W.gens[word]; init = W.one)

isLeftDescent(W, w, s) = s^w > W.data[:N]

function firstDescent(W, w)
    n, N = W.data[:rank], W.data[:N]
    return findfirst(s -> s^w > N, 1:n)
end

function coxeterWord(W, w)
    word = Int[]
    while !isidentity(w)
        a = firstDescent(W, w)
        push!(word, a)
        w = W.gens[a] * w
    end
    return word
end

reducedWord(W, word) = coxeterWord(W, permCoxeterWord(W, word))

function longestElt(W, J)
    wJ = W.one;  N = W.data[:N]
    while true
        i = findfirst(s -> s^wJ <= N, J)
        isnothing(i) && return wJ
        wJ = W.gens[J[i]] * wJ
    end
end

function prefixes(W, w)
    weak_quo(w, s) = isLeftDescent(W, w^-1, s) ? w * W.gens[s] : w
    return orbit(1:W.data[:rank], w, weak_quo)
end

function prefixes_with_edges(W, w)
    weak_quo(w, s) = isLeftDescent(W, w^-1, s) ? w * W.gens[s] : w
    orbit_with_edges(1:W.data[:rank], w, weak_quo)
end

longestCosetElt(W, J, L) = longestElt(W, J) * longestElt(W, L)

parabolicTransversal(W, J) = prefixes(W, longestCosetElt(W, J, 1:W.data[:rank]))

tackOn(x, s) = sort(union(x, [s]))

takeAway(x, s) = sort(setdiff(x, [s]))

onSortedTuples(tup, a) = sort([x^a for x in tup])

function shape(W, J)
    function onParabolics(K, s)
        return onSortedTuples(K, longestCosetElt(W, K, tackOn(K, s)))
    end
    sort(orbit(1:W.data[:rank], J, onParabolics))
end

function shapes(W)
    onShapes(x, s) = shape(W, takeAway(x[1], s))
    S = 1:W.data[:rank]
    orbit(S, shape(W, collect(S)), onShapes)
end

function shape_with_edges(W, J)
    function onParabolics(K, s)
        onSortedTuples(K, longestCosetElt(W, K, tackOn(K, s)))
    end
    orbit_with_edges(1:W.data[:rank], J, onParabolics)
end

function shape_with_transversal(W, J)
    S = 1:W.data[:rank]
    list = [J]
    index = Dict(J => 1)
    reps = [W.one]
    for (i, K) in enumerate(list)
        for s in setdiff(S, K)
            a = longestCosetElt(W, K, tackOn(K, s))
            L = onSortedTuples(K, a)
            get!(index, L) do
                push!(list, L);
                push!(reps, reps[i] * a);
                length(list)
            end
        end
    end
    return (list = list, reps = reps)
end

function parabolicComplement(W, J)
    S = 1:W.data[:rank]
    list = [J]
    index = Dict(J => 1)
    reps = [W.one]
    gens = (ears = Set(Perm[]), eyes = Set(Perm[]))
    for (i, K) in enumerate(list)
        for s in setdiff(S, K)
            a = longestCosetElt(W, K, tackOn(K, s))
            L = onSortedTuples(K, a)
            t = onRight(reps[i], a)
            l = get!(index, L) do
                push!(list, L)
                push!(reps, t)
                length(list)
            end
            if t != reps[l]
                push!(i == l ? gens.ears : gens.eyes, t / reps[l])
            end
        end
    end
    return gens
end

function minLenCons(W, x)
    list = [x]
    lx = coxeterLength(W, x)
    for (i, y) in enumerate(list)
        for s in W.gens
            z = y^s
            lz = coxeterLength(W, z)
            lz < lx && return minLenCons(W, z) # recurse!
            lz > lx && continue
            z in list || push!(list, z)
        end
    end
    return list
end

function coxeterMinRep(W, w)
    v = first(minLenCons(W, w))
    K = unique(coxeterWord(W, v))
    shape = shape_with_transversal(W, K)
    (L, j) = findmin(shape.list)
    return minimum(minLenCons(W, v^shape.reps[j]))
end

function coxeterConjugacyClasses(W)
    onMinReps(x, a) = coxeterMinRep(W, x * a)
    return orbit(reflections(W), W.one, onMinReps)
end

end # module
