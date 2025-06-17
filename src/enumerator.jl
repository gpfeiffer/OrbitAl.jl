#############################################################################
##
## enumerator.jl
##
module enumerator

using ..variants
using ..orbits
using ..permutation
using ..permgroup
import ..permgroup: PermGp

export Node, is_active, flat
export coset_table, compact_table

# ## Example

# * The complex reflection group $G_{12}$ has a presentation
# $$
# \langle
# s_1, s_2, s_3 \mid
# s_1^2 = s_2^2 = s_3^2 = 1,\,
# s_1 s_2 s_3 s_1 = s_2 s_3 s_1 s_2 = s_3 s_1 s_2 s_3
# \rangle
# $$
# * Let's try and enumerate its elements systematically.

# ## Smart Nodes
#
# * We will use a similar data structure, `Node`, for the purpose of coset enumeration.
# * Here, the `idx` attribute is used to identify `Node` objects.
# * And a `data` attribute is shared between all `Node` objects.

# %%
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

# %% [markdown]
# ## Coset Enumeration
#
# * Q: What is $G = \langle S \mid R \rangle$?
# * A: Todd-Coxeter!

# %% [markdown]
# * Suppose that a group $G$ is given by a **presentation** $\langle S \mid R \rangle$, consisting of a (finite) set $S$ of abstract **generators** $s_1, s_2, \dots, s_k$, and a (finite) list $R$ of **relations** $l_j = r_j$, where both $l_j$ and $r_j$ are words in $S \cup S^{-1}$.
#
# * For convenience, we assume that $S$ is closed under inverses: $S  = S^{-1}$.
#
# * We wish to enumerate the elements of $G$ (hoping that $G$ is a finite group), or more generally, the cosets of a subgroup $H$ of $G$ (hoping that $H$ has finite index in $G$).
#
# * A priori, neither the domain $X$ being acted upon (by $G$), nor the edges of the action graph are known.
#
# * Strategy: define new nodes as images of old nodes under a generator, but be prepared to identify this node with an existing one, if the relations imply they are the same.
#
# * For this, each `Node` object $x$ has
#   * a unique ID `idx` (where `idx` $ = n \iff x = x_n$),
#   * a word `word` $ \in S^*$ (corresponding to a path in the BFS spanning tree of the action graph),
#   * images $x$.`next`$[s] = x.s$ for each $s \in S$ (where $x.s \in X \cup \{ \perp \}$)
#   * a reference $x$.`flat`$ \in X \cup \{ \perp \}$ to the node it has possibly been replaced by.
#
# * Eventually, we want that x.s \in X for all x \in X, s \in S.

# * A node x is **active** if x.flat = \perp.
#
is_active(node::Node) = isnothing(node.flat)

# * Each node x \in X has an associated active node x^{\flat} defined recursively as x if x is active, and as (x.flat)^{\flat} otherwise
#
flat(node::Node) = is_active(node) ? node : flat(node.flat)

# * Recall that $S^{-1} = S$.  Assume that `data.invr` holds the map $s \mapsto s^{-1}$.
# * In words, we write $-s$ for $s^{-1}$.
# * So to find $x.s$ for $s \in S = S^{-1}$ we need to replace $s$ by `data.invr`$[-s]$ first, if $s < 0$.

# %%
function getImage(node::Node, s::Int)
    s < 0 ? node.next[node.data[:invr][-s]] : node.next[s]
end

# * To sprout a new node $x.s$:
#
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


# * Both actions need only be defined on the generators s \in S, and can then be applied to words in S.
#
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

# ### Enumerate!
#
# * We now formulate the `tabulate` procedure which takes a presentation `genrel` for a group $G$ as input and produces a permutation group as output.  Specifically, `genrel` has components
#   * `gens`: a list `[1..n]` of abstract generators $S = S^{-1}$
#   * `rels`: a list of relations expressed as pairs of word in $S$
#   * `invr`: the map $S \to S: s \mapsto s^{-1}$
#   * `sbgp`: a subset of $S$, generating a subgroup $H$ of $G$.
#
function trivialCoset(data, sbgp)
    node = Node([], data)
    for word in sbgp    # close the subgroup tables.
        trace(node, word)
    end
    return node
end


# ###  Tracing Words

# * To trace a node x under a word w means to make sure that x.w = x, using the sprouting action, i.e., creating new intermediate cosets as needed
# * If w \in H then x1.w = x1 should hold, for the trivial coset x1
# * If $l = r$ is a relation then $w:= l/r = 1$ and $x.w = x$ should hold for any $x \in X$.
# * In any case, before applying the last letter of w, we carefully check if the resulting coset is already known or not.
#
function trace(node::Node, word::Vector{Int})
    other = nodeUnderWordSprout(node, word[1:end-1])
    updateEdge(other, word[end], node)
end

# ### Processing a Node under a Generator
#
# * To find $x.s$, use variants of the relations to express $s$ as a word $w$ in the generators and check if $x.w$ is determined already.  If so, carefully set $x.s$ to $s.w$.  If this doesn't work out, create a new node $x.s$.
#
function process(node::Node, s::Int)
    for variant in node.data[:variants][s]
        if is_active(node)
            next = nodeUnderWordPartial(node, variant)
            isnothing(next) || updateEdge(node, s, next)
        end
    end
    if is_active(node) && isnothing(node.next[s])
        return sprout(node, s)
    else
        return flat(node).next[s]  # assuming that flat(node).next is done
    end
end

# ### Edges

# * In the (directed) graph of a group action, an edge $x \stackrel{s}{\longrightarrow} y$ always comes with the opposite edge $y \stackrel{s^{-1}}{\longrightarrow} x$.
# * Thus, carefully updating $x.s = y$ always refers to two edges of the graph.
#
function updateEdge(node::Node, s::Int, next::Node)
    setImage(node, s, next)
    setImage(next, node.data[:invr][s], node)
end

# * Carefully setting x.s to y means
#   * checking if x.s is already defined; if not, set x.s to y.
#   * Otherwise, with x.s = z, say, if y = z: there is nothing to do.
#   * Otherwise, merge z and y (keeping the older one active) and live with the consequenses ...
#
function setImage(node::Node, s::Int, next::Node)
    if isnothing(node.next[s])
        node.next[s] = next           # deduction!
    else
        y, z = flat(next), flat(node.next[s])
        y > z && ((y, z) = (z, y))        # sort
        y == z || mergeNodes(z, y)  # coincidence!
    end
end

# * to merge nodes z and y:
#   * set z.flat to y
#   * for each $z.s \neq {\perp}$, carefully update $z.s = y.s$.
#
function mergeNodes(node::Node, other::Node)
    node.flat = other
    node.data[:active] -= 1
    for (s, next) in enumerate(node.next)
        isnothing(next) || updateEdge(other, s, next)
    end
end

##  Finally
function coset_table(genrel, sbgp)
    data = Dict(
        :list => Node[],
        :active => 0,
        :gens => genrel.gens,
        :invr => genrel.invr
    )
    data[:variants] = variantsRelations(genrel)
    node = trivialCoset(data, sbgp)
    return orbitx(data[:gens], data[:list], process)
end

function compact_table(list)
    acti = filter(is_active, list)
    for (i, node) in enumerate(acti)
        node.next = flat.(node.next)
        node.idx = i
    end
    return acti
end

function PermGp(list::Vector{Node})
    imgs = [[next.idx for next in node.next] for node in compact_table(list)]
    gens = [Perm([x[i] for x in imgs]) for i in eachindex(imgs[1])]
    return PermGp(gens, one(gens[1]))
end

end # module
