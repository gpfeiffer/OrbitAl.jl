#############################################################################
##
#A  enumerator.jl                                                     OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  A simple, modular coset enumerator
##
module enumerator

using ..variants
using ..orbits
using ..permutation
using ..permgroup
import ..permgroup: PermGp

export Node, is_active, flat
export coset_table, compact_table

mutable struct Node
    idx::Int
    word::Vector{Int}
    flat::Union{Node, Nothing}
    next::Vector{Union{Node, Nothing}}
    data::Dict{Symbol, Any}
    function Node(word, data)
        l = length(data[:list])
        next = similar(data[:gens], Nothing)
        node = new(l + 1, word, nothing, next, data)
        push!(data[:list], node)
        data[:active] += 1
        return node
    end
end

## print a node
Base.show(io::IO, node::Node) = print(io, "Node(", node.idx, ")")
Base.show(io::IO, ::MIME"text/plain", node::Node) = print(io, "Node(", node.idx, ", word=", node.word, ")")

##  comparison
import Base: ==
==(node::Node, other::Node) = node.idx == other.idx
Base.hash(node::Node, h::UInt) = hash(node.idx, h)
Base.isless(node::Node, other::Node) = node.idx < other.idx


## A node x is **active** if x.flat = \perp.
is_active(node::Node) = isnothing(node.flat)

## Each node x \in X has an associated active node x^{\flat} defined
## recursively as x if x is active, and as (x.flat)^{\flat} otherwise
flat(node::Node) = is_active(node) ? node : flat(node.flat)

##  get node.next[s], allowing for inverses
function getImage(node::Node, s::Int)
    s < 0 ? node.next[node.data[:invr][-s]] : node.next[s]
end

## To sprout a new node $x.s$:
function sprout(node::Node, s::Int)
    @assert node.next[s] === nothing "node.$s already defined"
    next = Node(onWords(node.word, s), node.data) # new node
    node.next[s] = next
    next.next[node.data[:invr][s]] = node
    return next
end

# We will work with two distinct actions:
# * a **partial action** which returns `nothing` for undefined images
# * a **sprouting action** which sprouts a new node if necessary.
function onNodes(node::Node, s::Int, func::Function)
    next = getImage(node, s)
    isnothing(next) ? func(node, s) : flat(next)
end
onNodesPartial(node, s) = onNodes(node, s, (x, a) -> nothing)
onNodesSprout(node, s) = onNodes(node, s, sprout)

# Extend both actions to words in S
function nodeUnderWordSprout(node::Node, word::Vector{Int})
    for s in word
        node = onNodesSprout(node, s)
    end
    return node
end

function nodeUnderWordPartial(node::Node, word::Vector{Int})
    for s in word
        node = onNodesPartial(node, s)
        isnothing(node) && return node
    end
    return node
end

## construct the trivial coset and close subgroup tables.
function trivialCoset(data, sbgp)
    node = Node([], data)
    for word in sbgp    # close the subgroup tables.
        trace(node, word)
    end
    return node
end

# Tracing node x under word w ensures x.w = x, using the sprouting action
function trace(node::Node, word::Vector{Int})
    other = nodeUnderWordSprout(node, word[1:end-1])
    updateEdge(other, word[end], node)
end

# To find x.s, try all variants of the relations, create x.s if that fails.
function finalize(node::Node, s::Int)
    for variant in node.data[:variants][s]
        if is_active(node)
            next = nodeUnderWordPartial(node, variant)
            isnothing(next) || updateEdge(node, s, next)
        end
    end
    if is_active(node) && isnothing(node.next[s])
        return sprout(node, s)
    end
    return flat(node).next[s]  # assuming that flat(node).next is done
end

# carefully update the edge x.s = y in both directions
function updateEdge(node::Node, s::Int, next::Node)
    setImage(node, s, next)
    setImage(next, node.data[:invr][s], node)
end

# Carefully set x.s to y making deductions and stacking coincidences
function setImage(node::Node, s::Int, next::Node)
    if isnothing(node.next[s])
        node.next[s] = next           # deduction!
    else
        y, z = flat(next), flat(node.next[s])
        y > z && ((y, z) = (z, y))        # sort
        y == z || mergeNodes(z, y)  # coincidence!
    end
end

# merge nodes z and y, keeping the older one
function mergeNodes(node::Node, other::Node)
    node.flat = other
    node.data[:active] -= 1
    for (s, next) in enumerate(node.next)
        isnothing(next) || updateEdge(other, s, next)
    end
end

## construct the coset table as an orbit
function coset_table(genrel, sbgp)
    data = Dict(
        :list => Node[],
        :active => 0,
        :gens => genrel.gens,
        :invr => genrel.invr
    )
    data[:variants] = variantsRelations(genrel)
    node = trivialCoset(data, sbgp)
    return orbitx(data[:gens], data[:list], finalize)
end

# drop redundant nodes and relabel.
function compact_table(list)
    acti = filter(is_active, list)
    for (i, node) in enumerate(acti)
        node.next = flat.(node.next)
        node.idx = i
    end
    return acti
end

# how to convert coset table into a perm group
function PermGp(list::Vector{Node})
    imgs = [[next.idx for next in node.next] for node in compact_table(list)]
    gens = [Perm([x[i] for x in imgs]) for i in eachindex(imgs[1])]
    return PermGp(gens, one(gens[1]))
end

end # module
