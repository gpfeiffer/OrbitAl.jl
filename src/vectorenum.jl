#############################################################################
##
#A  vectorenum.jl                                                     OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  A vector enumerator for Hecke algebras: the coset enumerator, with
#C  cosets as basis vectors and linear Union-Find for the coincidences
##
module vectorenum

using ..sparsevec
import ..unionfind: find, unite!, isunit
import ..enumerator: is_active, active_cosets

export VectorTable, vector_enum, knownact, matrices, is_active, active_cosets

const Row{C} = Union{Nothing, SparseVec{C}}

"""
    VectorTable{C}

A vector table under construction, with coefficients of type `C`.  The cosets
`1, 2, 3, ...` are basis vectors ``e_1, e_2, e_3, \\dots``; `next[x][s]` is the
image ``e_x.s``, as a `SparseVec`, or `nothing` if not yet known.  `parent`
holds the linear relations between the cosets (see
[`OrbitAl.unionfind`](@ref Union-Find)), `dead` the rows `(i, s, w)` of
eliminated cosets, which still have to be satisfied, `pending` the relations
that have no unit coefficient (yet), and `words[x]` is the word that defined
coset `x`.
"""
mutable struct VectorTable{C}
    one::C                                     # the coefficient 1
    next::Vector{Vector{Row{C}}}               # next[x][s]: x.s, or nothing
    words::Vector{Vector{Int}}                 # the word that defined each coset
    parent::Dict{Int, SparseVec{C}}            # linear Union-Find on the cosets
    dead::Vector{Tuple{Int, Int, SparseVec{C}}}    # rows (i, s, w) of eliminated cosets
    pending::Vector{SparseVec{C}}              # relations without a unit coefficient
    relators::Vector{Vector{Tuple{C, Vector{Int}}}}
    ngens::Int
end

VectorTable(c::C, relators, ngens) where C =
    VectorTable{C}(c, [], [], Dict(), [], [], relators, ngens)

basis(T::VectorTable, i) = SparseVec([i], [T.one])

"""
    is_active(T::VectorTable, x)

Whether coset `x` is active, that is, has not been eliminated by a relation.
"""
is_active(T::VectorTable, x) = !haskey(T.parent, x)

"""
    active_cosets(T::VectorTable)

The list of active cosets of `T`.  Once the enumeration is complete, these
form a basis of the module.
"""
active_cosets(T::VectorTable) = filter(x -> is_active(T, x), eachindex(T.next))

function newcoset!(T::VectorTable, word)
    push!(T.next, fill(nothing, T.ngens))
    push!(T.words, word)
    return length(T.next)
end

##  the quadratic relators T_s^2 - (q-1) T_s - q, and l - r for the relations
##  l = r of the presentation, as lists of pairs (coefficient, word)
function relators(genrel, q, c)
    rels = [[(c, [s, s]), (-(q - 1) * c, [s]), (-q * c, Int[])] for s in genrel.gens]
    append!(rels, [[(c, l), (-c, r)] for (l, r) in genrel.rels])
end

"""
    knownact(T::VectorTable, v, word)

The image `v.word`, if the table knows it, and `nothing` otherwise.
"""
function knownact(T::VectorTable, v, word)
    for s in word
        v = find(T.parent, v)
        r = zero(v)
        for (j, h) in zip(v.poss, v.vals)
            isnothing(T.next[j][s]) && return nothing
            r += h * T.next[j][s]
        end
        v = r
    end
    return find(T.parent, v)
end

##  terms (u, s), standing for u.s, plus a constant vector: a relation,
##  a deduction, or nothing yet
function evaluate!(T::VectorTable, terms, constant)
    v, unknown = constant, Dict{Tuple{Int, Int}, typeof(T.one)}()
    for (u, s) in terms
        u = find(T.parent, u)
        for (j, h) in zip(u.poss, u.vals)
            w = T.next[j][s]
            isnothing(w) ? (unknown[(j, s)] = get(unknown, (j, s), zero(h)) + h) : (v += h * w)
        end
    end
    filter!(p -> !iszero(p.second), unknown)
    v = find(T.parent, v)
    if isempty(unknown)
        return iszero(v) ? :none : relation!(T, v)
    elseif length(unknown) == 1 && isunit(only(values(unknown)))
        ((j, s), h) = only(unknown)
        T.next[j][s] = find(T.parent, -(inv(h) * v))         # a deduction
        return :deduced
    end
    return :wait
end

##  a relation v = 0 eliminates a coset, whose row becomes a list of equations,
##  or waits for a unit coefficient
function relation!(T::VectorTable, v)
    i = unite!(T.parent, v, zero(v))
    i == -1 && (push!(T.pending, v); return :wait)
    i > 0 || return :none
    for s in 1:T.ngens
        w = T.next[i][s]
        isnothing(w) || push!(T.dead, (i, s, w))
        T.next[i][s] = nothing
    end
    return :relation
end

##  a relator at coset x, prepared for evaluate!, or nothing if an entry
##  before a last letter is unknown
function terms_at(T::VectorTable, x, rel)
    terms, constant = Tuple{SparseVec{typeof(T.one)}, Int}[], zero(basis(T, x))
    for (c, word) in rel
        if isempty(word)
            constant += c * basis(T, x)
        else
            u = knownact(T, basis(T, x), word[1:end-1])
            isnothing(u) && return nothing
            push!(terms, (c * u, word[end]))
        end
    end
    return terms, constant
end

##  apply all relators at all active cosets, the rows of eliminated cosets,
##  and the pending relations, until nothing changes
function close!(T::VectorTable)
    changed = true
    while changed
        changed = false
        for x in active_cosets(T), rel in T.relators
            is_active(T, x) || continue
            tc = terms_at(T, x, rel)
            isnothing(tc) && continue
            evaluate!(T, tc...) in (:relation, :deduced) && (changed = true)
        end
        for k in reverse(eachindex(T.dead))
            (i, s, w) = T.dead[k]
            evaluate!(T, [(T.parent[i], s)], -w) == :wait || (deleteat!(T.dead, k); changed = true)
        end
        pending, T.pending = T.pending, empty(T.pending)
        for v in pending                       # retry, with what is known now
            relation!(T, v) == :relation && (changed = true)
        end
    end
end

"""
    vector_enum(genrel, q, J, Jcoef; defs = [])

Enumerate a basis of the module ``H \\otimes_{H_J} M`` of the Hecke algebra
``H`` with parameter `q` of the presentation `genrel` (with fields `gens` and
`rels`, the relations `l = r` holding in ``H`` for the ``T_s``, alongside the
quadratic relations ``T_s^2 = (q-1)T_s + q``).  The generators `J` act on the
first coset ``x_1`` as the coefficients `Jcoef`: scalars, such as `q`, or
elements of the Hecke algebra ``H_J``, such as [`HElt`](@ref
OrbitAl.hecke.HElt)s, acting from the left.

New cosets are defined in breadth-first order, or along the words `defs`
first, where possible.  Returns a [`VectorTable`](@ref
OrbitAl.vectorenum.VectorTable).  A relation without a unit coefficient
waits; if any are left at the end, the active cosets span the module but need
not be a basis.
"""
function vector_enum(genrel, q, J, Jcoef; defs = Vector{Int}[])
    c = isempty(Jcoef) ? one(q) : one(first(Jcoef))
    T = VectorTable(c, relators(genrel, q, c), length(genrel.gens))
    x = newcoset!(T, Int[])
    for (t, h) in zip(J, Jcoef)
        T.next[x][t] = h * basis(T, x)
    end
    while true
        close!(T)
        open = [(x, s) for x in active_cosets(T) for s in 1:T.ngens if isnothing(T.next[x][s])]
        isempty(open) && break
        for w in defs                          # a prescribed definition, if possible
            k = findfirst(==(w[1:end-1]), T.words)
            if !isnothing(k) && (k, w[end]) in open
                open = [(k, w[end])]; break
            end
        end
        (x, s) = first(open)
        T.next[x][s] = basis(T, newcoset!(T, [T.words[x]; s]))
    end
    isempty(T.pending) || @warn "$(length(T.pending)) relations without a unit coefficient left"
    return T
end

"""
    matrices(T::VectorTable)

The matrices of the generators on the active cosets, renumbered `1, ..., n`:
row `k` of the matrix of `s` holds the coefficients of ``x_k.s``.
"""
function matrices(T::VectorTable)
    acti = active_cosets(T)
    idx = Dict(x => k for (k, x) in enumerate(acti))
    mats = [fill(zero(T.one), length(acti), length(acti)) for _ in 1:T.ngens]
    for (k, x) in enumerate(acti), s in 1:T.ngens
        v = find(T.parent, T.next[x][s])
        for (j, h) in zip(v.poss, v.vals)
            mats[s][k, idx[j]] = h
        end
    end
    return mats
end

function Base.show(io::IO, ::MIME"text/plain", T::VectorTable)
    acti = active_cosets(T)
    num = Dict(x => k for (k, x) in enumerate(acti))
    println(io, "VectorTable: ", length(acti), " active of ", length(T.next), " cosets")
    cf(h) = isone(h) ? "" : "(" * sprint(show, h; context = :typeinfo => typeof(h)) * ")"
    function vstr(v)
        isnothing(v) && return "?"
        v = find(T.parent, v)
        iszero(v) && return "0"
        join([cf(h) * "x" * string(num[j]) for (j, h) in zip(v.poss, v.vals)], " + ")
    end
    for x in acti[1:min(end, 8)]
        println(io, "x", num[x], isempty(T.words[x]) ? "" : " = x1." * join(T.words[x], "."))
        for s in 1:T.ngens
            println(io, "    .", s, " = ", vstr(T.next[x][s]))
        end
    end
    length(acti) > 8 && print(io, "  ...")
end

end # module
