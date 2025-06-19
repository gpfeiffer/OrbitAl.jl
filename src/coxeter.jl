#############################################################################
##
module coxeter

# using LinearAlgebra

using ..permutation
using ..orbits
using ..permgroup

import Base: ^
import ..permutation: Perm
import ..permutation: isidentity, last_moved
import ..permgroup: PermGp

export coxeterGraph, cartanMat, CoxeterGp, coxeterConjugacyClasses
export coxeterLength, coxeterWord, reflections
export shapes

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

function orbits_with_words_and_edges(aaa, xxx, under)
    words = Array{Int}[[i] for i in eachindex(xxx)]
    list = copy(xxx)
    edges = []
    for (i, y) in enumerate(list)
        for (k, a) in enumerate(aaa)
            z = under(y, a)
            l = findfirst(==(z), list)
            if isnothing(l)
                push!(list, z)
                push!(words, onWords(words[i], k))
                l = length(list)
            end
            push!(edges, (i, l))
        end
    end
    return (list = list, edges = edges, words = words)
end

Perm(a, xxx, under) = Perm(indexin([under(x, a) for x in xxx], xxx))

struct CoxeterGp
    gens::Vector{Perm}
    one::Perm
    data::Dict{Symbol, Any}
#    CoxeterGp(gens, one) = new(gens, one, Dict())
end

function reflection(C::Matrix{Int}, s::Int)
    mat = C^0
    mat[:,s] = mat[s,:] - C[s,:]
    return mat
end

function CoxeterGp(C)
    m1 = C^0
    S = axes(C,1)
    mats = [reflection(C, s) for s in S]
    roots = orbitx_with_words(mats, [m1[i:i,:] for i in S], onRoots)
    data = Dict(:mats => mats, :roots => roots, :rank => length(S))
    data[:N] = length(roots.list)
    data[:phi] = [roots.list; -roots.list]
    data[:perms] = [Perm(m, data[:phi], onRight) for m in mats]
    return CoxeterGp(data[:perms], data[:perms][1]^0, data)
end

PermGp(group::CoxeterGp) = PermGp(group.gens, group.one)

size(group::CoxeterGp) = sizeOfGroup(PermGp(group))

^(a::Perm, group::CoxeterGp) = a^PermGp(group)

data(group::CoxeterGp) = group.data

function reflections(W::CoxeterGp)
    refl(w) = W.gens[w[1]]^prod(W.gens[w[2:end]]; init=W.one)
    refl.(data(W)[:roots].words)
end

function coxeterLength(W, w)
    N = data(W)[:N]
    count(i^w > N for i in 1:N)
end

permCoxeterWord(W, word) = prod(W.gens[word]; init = W.one)

isLeftDescent(W, w, s) = s^w > data(W)[:N]

function firstDescent(W, w)
    n, N = data(W)[:rank], data(W)[:N]
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

function longestElement(W, J)
    wJ = W.one
    N = data(W)[:N]
    J = collect(J)
    while true
        i = findfirst(s -> s^wJ <= N, J)
        isnothing(i) && return wJ
        wJ = W.gens[J[i]] * wJ
    end
end

function prefixes(W, w)
    onRightDescents(w, s) = isLeftDescent(W, w^-1, s) ? w * W.gens[s] : w
    return orbit(1:data(W)[:rank], w, onRightDescents)
end

function prefixes_with_edges(W, w)
    onRightDescents(w, s) = isLeftDescent(W, w^-1, s) ? w * W.gens[s] : w
    orbit_with_edges(1:data(W)[:rank], w, onRightDescents)
end

longestCosetElement(W, J, L) = longestElement(W, J) * longestElement(W, L)

parabolicTransversal(W, J) = prefixes(W, longestCosetElement(W, J, 1:data(W)[:rank]))

tackOn(x, s) = sort(union(x, s))

takeAway(x, s) = sort(setdiff(x, s))

onSortedTuples(tup, a) = sort([x^a for x in tup])

function shape(W, J)
    function onParabolics(K, s)
        return onSortedTuples(K, longestCosetElement(W, K, tackOn(K, s)))
    end
    sort(orbit(1:data(W)[:rank], J, onParabolics))
end

function shapes(W)
    onShapes(x, s) = shape(W, takeAway(x[1], s))
    S = 1:data(W)[:rank]
    orbit(S, shape(W, collect(S)), onShapes)
end

function shape_with_edges(W, J)
    function onParabolics(K, s)
        onSortedTuples(K, longestCosetElement(W, K, tackOn(K, s)))
    end
    orbit_with_edges(1:data(W)[:rank], J, onParabolics)
end

function shape_with_transversal(W, J)
    S = 1:data(W)[:rank]
    list = [J]
    index = Dict(J => 1)
    reps = [W.one]
    for (i, K) in enumerate(list)
        for s in setdiff(S, K)
            a = longestCosetElement(W, K, tackOn(K, s))
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
    S = 1:data(W)[:rank]
    list = [J]
    index = Dict(J => 1)
    reps = [W.one]
    gens = (ears = Set(Perm[]), eyes = Set(Perm[]))
    for (i, K) in enumerate(list)
        for s in setdiff(S, K)
            a = longestCosetElement(W, K, tackOn(K, s))
            L = onSortedTuples(K, a)
            t = onRight(reps[i], a)
            j = get!(index, L) do
                push!(list, L)
                push!(reps, t)
                length(list)
            end
            if t != reps[j]
                push!(i == j ? gens.ears : gens.eyes, t /reps[j])
            end
        end
    end
    return gens
end

function minConjugates(W, x)
    list = [x]
    lx = coxeterLength(W, x)
    i = 0
    while i < length(list)
        i += 1;
        y = list[i]
        for s in W.gens
            z = y^s
            lz = coxeterLength(W, z)
            if lz == lx
                z in list || push!(list, z)
            elseif lz < lx  # reset list
                list = [z]
                lx = lz
                i = 0
                break
            end
        end
    end
    return list
end

is_trivial(group::CoxeterGp) = all(isidentity, group.gens)
last_moved(group::CoxeterGp) = max(last_moved.(group.gens)...)

function coxeterMinRep(W, w)
    v = first(minConjugates(W, w))
    K = unique(coxeterWord(W, v))
    shape = shape_with_transversal(W, K)
    (L, j) = findmin(shape.list)
    return minimum(minConjugates(W, v^shape.reps[j]))
end

function coxeterConjugacyClasses(W)
    onMinReps(x, a) = coxeterMinRep(W, x * a)
    return orbit(reflections(W), W.one, onMinReps)
end

end # module
