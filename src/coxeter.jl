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

export coxeterGraph, coxeterMat, cartanMat, CoxeterGp, coxeterConjugacyClasses
export coxeterLength, coxeterWord, permCoxeterWord, reflections
export size, prefixes, prefixes_with_edges, shapes
export longestElt, parabolicTransversal

"""
    coxeterGraph(series, rank)

Return the edge list of the Coxeter diagram for the group of given `series`
("A", "B", "C", "D", "E", "F", ...) and `rank`.
Each edge `(i, j)` indicates that generators `s_i` and `s_j` have Coxeter order 3.

# Examples
```jldoctest
julia> using OrbitAl

julia> coxeterGraph("A", 3)
2-element Vector{Tuple{Int64, Int64}}:
 (1, 2)
 (2, 3)
```
"""
function coxeterGraph(series::String, rank::Int)
    edges = [(j-1,j) for j in 2:rank]     # type A: chain
    series >= "D" && (edges[1] = (1,3))   # type D: add fork
    series >= "E" && (edges[2] = (2,4))   # type E: another fork
    return edges
end

"""
    coxeterMat(series, rank)

Return the Coxeter matrix for the group of given `series` and `rank`.
Entry `m[i,j]` is the order of the product `s_i * s_j` of generators.

# Examples
```jldoctest
julia> using OrbitAl

julia> coxeterMat("A", 2)
2×2 Matrix{Int64}:
 1  3
 3  1
```
"""
function coxeterMat(series::String, n::Int)
    mat = zeros(Int, n, n)^0 # identity mat
    edges = coxeterGraph(series, n)
    for j in 1:n
        for i in 1:j-1
            mat[i,j] = mat[j,i] = (i,j) in edges ? 3 : 2
        end
    end
    return mat
end


mij(i, j, m) = repeat([i,j], div(m+1,2))[1:m]

function coxeterPresentation(series::String, rank::Int)
    gens = collect(1:rank)
    genrel = (gens = gens, invr = gens, rels = [])
    edges = coxeterGraph(series, rank)
    for j in gens
        for i in 1:j-1
            if (i, j) in edges
                push!(genrel.rels, [[i,j,i],[j,i,j]])
            else
                push!(genrel.rels, [[i,j],[j,i]])
            end
        end
    end
    return genrel
end

"""
    cartanMat(series, rank)

Return the Cartan matrix for the Lie algebra of given `series` and `rank`.

# Examples
```jldoctest
julia> using OrbitAl

julia> cartanMat("A", 2)
2×2 Matrix{Int64}:
  2  -1
 -1   2
```
"""
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

"""
    CoxeterGp(C)

Construct a Coxeter group from the Cartan matrix `C`, realized as a permutation
group acting on its root system. The generators correspond to simple reflections.

# Examples
```jldoctest
julia> using OrbitAl

julia> W = CoxeterGp(cartanMat("A", 2));

julia> sizeOfGroup(W)
6
```
"""
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

"""
    reflections(W)

Return all reflections in the Coxeter group `W` as a vector of permutations.
Reflections are conjugates of the generators; their count equals the number of
positive roots.

# Examples
```jldoctest
julia> using OrbitAl

julia> W = CoxeterGp(cartanMat("A", 2));

julia> length(reflections(W))
3
```
"""
function reflections(W)
    refl(w) = W.gens[w[1]]^prod(W.gens[w[2:end]]; init=W.one)
    refl.(W.data[:roots].words)
end

"""
    coxeterLength(W, w)

Return the Coxeter length of `w` in `W`, i.e., the number of positive roots
sent to negative roots by `w`, which equals the length of a shortest
generator word for `w`.

# Examples
```jldoctest
julia> using OrbitAl

julia> W = CoxeterGp(cartanMat("A", 2));

julia> coxeterLength(W, W.gens[1] * W.gens[2])
2
```
"""
function coxeterLength(W, w)
    N = W.data[:N]
    return count(i^w > N for i in 1:N)
end

"""
    permCoxeterWord(W, word)

Return the permutation in `W` corresponding to the word `word`,
a vector of generator indices.
"""
permCoxeterWord(W, word) = prod(W.gens[word]; init = W.one)

isLeftDescent(W, w, s) = s^w > W.data[:N]

function firstDescent(W, w)
    n, N = W.data[:rank], W.data[:N]
    return findfirst(s -> s^w > N, 1:n)
end

"""
    coxeterWord(W, w)

Return a reduced word for `w` in `W` as a vector of generator indices,
computed by iteratively stripping left descents.

# Examples
```jldoctest
julia> using OrbitAl

julia> W = CoxeterGp(cartanMat("A", 2));

julia> length(coxeterWord(W, W.gens[1] * W.gens[2]))
2
```
"""
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

"""
    longestElt(W, J)

Return the longest element of the parabolic subgroup of `W` generated by
the simple reflections indexed by `J`.
"""
function longestElt(W, J)
    wJ = W.one;  N = W.data[:N]
    while true
        i = findfirst(s -> s^wJ <= N, J)
        isnothing(i) && return wJ
        wJ = W.gens[J[i]] * wJ
    end
end

"""
    prefixes(W, w)

Return all right-sided prefixes of `w` in the right weak order of `W`,
i.e., the elements obtained from `w` by removing right descents one at a time.
"""
function prefixes(W, w)
    weak_quo(w, s) = isLeftDescent(W, w^-1, s) ? w * W.gens[s] : w
    return orbit(1:W.data[:rank], w, weak_quo)
end

"""
    prefixes_with_edges(W, w)

Like `prefixes`, but also returns the labeled edges of the Hasse diagram of
the right weak order interval below `w`.
"""
function prefixes_with_edges(W, w)
    weak_quo(w, s) = isLeftDescent(W, w^-1, s) ? w * W.gens[s] : w
    orbit_with_edges(1:W.data[:rank], w, weak_quo)
end

longestCosetElt(W, J, L) = longestElt(W, J) * longestElt(W, L)

"""
    parabolicTransversal(W, J)

Return the shortest left coset representatives for the parabolic subgroup of
`W` generated by `J`.
"""
parabolicTransversal(W, J) = prefixes(W, longestCosetElt(W, J, 1:W.data[:rank]))

tackOn(x, s) = sort(union(x, [s]))

takeAway(x, s) = sort(setdiff(x, [s]))

onSortedTuples(tup, a) = sort([x^a for x in tup])

"""
    shape(W, J)

Return the orbit of the parabolic subset `J` under conjugation by double coset
representatives, i.e., the set of all parabolic subsets conjugate to `J` in `W`.
"""
function shape(W, J)
    function onParabolics(K, s)
        return onSortedTuples(K, longestCosetElt(W, K, tackOn(K, s)))
    end
    sort(orbit(1:W.data[:rank], J, onParabolics))
end

"""
    shapes(W)

Return all parabolic shapes in `W`, i.e., the orbits of all parabolic subsets
under conjugation by double coset representatives.
"""
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

"""
    coxeterConjugacyClasses(W)

Return a list of minimum-length conjugacy class representatives in `W`,
one per conjugacy class.
"""
function coxeterConjugacyClasses(W)
    onMinReps(x, a) = coxeterMinRep(W, x * a)
    return orbitl(reflections(W), W.one, onMinReps)
end

end # module
