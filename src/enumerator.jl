#############################################################################
##
#A  enumerator.jl                                                     OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  A simple coset enumerator: orbit algorithm + Union-Find + relation variants,
#C  and the Reidemeister-Schreier method
##
module enumerator

using ..variants
using ..permutation
import ..orbits: orbit_with_words, onPoints
import ..permgroup: PermGp
import ..unionfind: find, unite!

export CosetTable, coset_table, is_active, active_cosets, perms
export schreier_edges, rewrite, reidemeister_schreier

"""
    CosetTable

A coset table under construction.  Cosets are numbered `1, 2, 3, ...`;
`next[x][s]` is the image of coset `x` under generator `s` (`0` if not yet
known), `parent` is a Union-Find forest on the cosets (`parent[x] == x` if `x`
is active), `queue` holds the pending coincidences, and `words[x]` is the word that defined coset `x`.
"""
mutable struct CosetTable
    next::Vector{Vector{Int}}             # next[x][s], 0 = not yet known
    parent::Vector{Int}                   # the Union-Find forest on the cosets
    queue::Vector{Tuple{Int,Int}}         # the pending coincidences
    words::Vector{Vector{Int}}            # the word that defined each coset
    invr::Vector{Int}                     # invr[s] is the inverse of s
    variants::Vector{Vector{Vector{Int}}} # the relation variants
    active::Int                           # the number of active cosets
end

CosetTable(genrel) = CosetTable([], [], [], [], genrel.invr, variantsRelations(genrel), 0)

##  Union-Find: the active coset that x has been merged into
find(T::CosetTable, x) = find(T.parent, x)

"""
    is_active(T::CosetTable, x)

Whether coset `x` is active, that is, has not been merged into another coset.
"""
is_active(T::CosetTable, x) = T.parent[x] == x

"""
    active_cosets(T::CosetTable)

The list of active cosets of `T`.  Once the enumeration is complete, these
are the cosets of the subgroup, and their number is `T.active`.
"""
active_cosets(T::CosetTable) = filter(x -> is_active(T, x), eachindex(T.parent))

##  x.s as an active coset, or 0 if not yet known
image(T::CosetTable, x, s) = (y = T.next[x][s]; y > 0 ? find(T, y) : 0)

function newcoset!(T::CosetTable, word)
    push!(T.next, zeros(Int, length(T.invr)))
    push!(T.parent, length(T.parent) + 1)
    push!(T.words, word)
    T.active += 1
    return length(T.parent)
end

##  define a new coset y = x.s
function sprout!(T::CosetTable, x, s)
    y = newcoset!(T, [T.words[x]; s])
    T.next[x][s] = y
    T.next[y][T.invr[s]] = x
    return y
end

function Base.show(io::IO, ::MIME"text/plain", T::CosetTable)
    acti = active_cosets(T)
    println(io, "CosetTable: ", length(acti), " active of ", length(T.parent), " cosets")
    for x in acti[1:min(end, 12)]
        println(io, lpad(x, 4), ": ", rpad(join(T.words[x]), 10), [image(T, x, s) for s in eachindex(T.invr)])
    end
    length(acti) > 12 && print(io, "   ...")
end

##  the partial action gives up (returns 0) at an unknown image
function underWordPartial(T::CosetTable, x, word)
    for s in word
        x = image(T, x, s)
        x > 0 || return 0
    end
    return x
end

##  the defining action defines new cosets where necessary
function underWordSprout!(T::CosetTable, x, word)
    for s in word
        y = image(T, x, s)
        x = y > 0 ? y : sprout!(T, x, s)
    end
    return x
end

##  x.s = y, or a coincidence if x.s is already set
function setImage!(T::CosetTable, x, s, y)
    x, y = find(T, x), find(T, y)
    z = T.next[x][s]
    z == 0 ? (T.next[x][s] = y) : push!(T.queue, (y, z))
end

##  x.s = y and y.s^-1 = x
function link!(T::CosetTable, x, s, y)
    setImage!(T, x, s, y)
    setImage!(T, y, T.invr[s], x)
end

##  merge the cosets a and b, and their rows: Union-Find with consequences
function unite!(T::CosetTable, a, b)
    b = unite!(T.parent, a, b)   # the coset that is no longer active, or 0
    b > 0 || return b
    a = T.parent[b]              # the coset that b is merged into
    T.active -= 1
    for (t, y) in enumerate(T.next[b])     # move b's edges over to a
        y == 0 || link!(T, a, t, y)
    end
    return b
end

##  set x.s = y, and process the resulting coincidences
function updateEdge!(T::CosetTable, x, s, y)
    link!(T, x, s, y)
    for (a, b) in T.queue                      # the queue grows while we loop
        unite!(T, a, b)
    end
    empty!(T.queue)                            # ready for the next edge
end

##  make sure that x.word = x
function trace!(T::CosetTable, x, word)
    y = underWordSprout!(T, x, word[1:end-1])
    updateEdge!(T, y, word[end], x)
end

##  find x.s from the relation variants, or define it
function finalize!(T::CosetTable, x, s)
    for variant in T.variants[s]
        is_active(T, x) || return
        y = underWordPartial(T, x, variant)
        y == 0 || updateEdge!(T, x, s, y)
    end
    is_active(T, x) || return
    T.next[x][s] > 0 || sprout!(T, x, s)
end

"""
    coset_table(genrel, sbgp)

Enumerate the cosets of the subgroup generated by the words `sbgp` in the
group with presentation `genrel` (with fields `gens`, `invr` and `rels`).
Returns a `CosetTable`; its active cosets are the cosets of the subgroup.
"""
function coset_table(genrel, sbgp)
    T = CosetTable(genrel)
    x = newcoset!(T, Int[])
    for word in sbgp
        trace!(T, x, word)
    end
    x = 0
    while x < length(T.parent)                # the table grows while we loop
        x += 1
        for s in genrel.gens
            finalize!(T, x, s)
        end
    end
    return T
end

"""
    perms(T::CosetTable)

The permutations of the generators on the active cosets, renumbered `1, ..., n`.
"""
function perms(T::CosetTable)
    acti = active_cosets(T)
    idx = Dict(x => i for (i, x) in enumerate(acti))
    [Perm([idx[image(T, x, s)] for x in acti]) for s in eachindex(T.invr)]
end

##  how to convert a coset table into a permutation group
function PermGp(T::CosetTable)
    gens = perms(T)
    return PermGp(gens, one(gens[1]))
end

#############################################################################
##
##  Presentations of subgroups: the Reidemeister-Schreier method
##

"""
    schreier_edges(genrel, P, list)

The Schreier generators of the subgroup: the redundant pairs of Union-Find on
the coset graph, where `P` holds the permutations of the generators on the
cosets and `list` holds the cosets in BFS order, so that the forest is the
BFS tree.  Returns `label, edges, rels`: `label[(x, s)]` is `k` or `-k` for
the pair of non-tree edges `edges[k]`, and `0` for a tree edge; `rels` holds
a relator `[k, k]` for each `edges[k]` that is a loop of an involution.
"""
function schreier_edges(genrel, P, list)
    invr = genrel.invr
    parent = collect(eachindex(list))
    label = Dict{Tuple{Int,Int},Int}()
    edges = Tuple{Int,Int}[]
    rels = Vector{Vector{Int}}()
    for x in list, s in genrel.gens
        haskey(label, (x, s)) && continue
        y = x^P[s]
        if unite!(parent, x, y) > 0         # a tree edge, and its partner
            label[(x, s)] = label[(y, invr[s])] = 0
            continue
        end
        push!(edges, (x, s)); k = length(edges)
        label[(x, s)] = k
        if (y, invr[s]) == (x, s)
            push!(rels, [k, k])       # an involution
        else
            label[(y, invr[s])] = -k
        end
    end
    return label, edges, rels
end

"""
    rewrite(P, label, x, w)

The word `w`, read from the coset `x`, rewritten in the Schreier generators
labelled by `label` (see `schreier_edges`): the word for `t_x w t_{x.w}^-1`.
"""
function rewrite(P, label, x, w)
    rel = Int[]
    for s in w
        k = label[(x, s)]
        k == 0 || push!(rel, k)
        x = x^P[s]
    end
    return rel
end

"""
    reidemeister_schreier(genrel, sbgp)

A presentation of the subgroup generated by the words `sbgp` in the group
with presentation `genrel`, which must have finite index.  Returns a named
tuple: `index`, the Schreier transversal `word` (a word for each coset), the
generators `edges` (a non-tree edge `(x, s)` for the Schreier generator
`t_x s t_{x.s}^-1`), the relators `rels` (words in `±k` for `edges[k]`), and
the permutations `perms` of the generators on the cosets.
"""
function reidemeister_schreier(genrel, sbgp)
    P = perms(coset_table(genrel, sbgp))
    tvl = orbit_with_words(P[genrel.gens], 1, onPoints)     # a Schreier transversal
    label, edges, rels = schreier_edges(genrel, P, tvl.list)
    relators = [[l; reverse(genrel.invr[r])] for (l, r) in genrel.rels]
    for x in tvl.list, w in relators
        rel = rewrite(P, label, x, w)
        isempty(rel) || push!(rels, rel)
    end
    return (index = length(tvl.list), word = Dict(zip(tvl.list, tvl.words)), edges = edges, rels = rels, perms = P)
end

end # module
