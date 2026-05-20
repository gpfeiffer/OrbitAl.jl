#############################################################################
##
#A  simsgroup.jl                                                      OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  A simple minded implementation of the Schreier-Sims algorithm
##
module simsgroup

import Base: in, size, rand, intersect

using ..orbits
import ..permutation: Perm, isidentity, last_moved

export SimsGp, orbit_sims, cube

"""
    SimsGp(gens, one)

A permutation group represented by a Schreier-Sims stabilizer chain.
`gens` is a list of generating permutations; `one` is the group identity.
The stabilizer chain is built lazily on first use.
"""
mutable struct SimsGp
    gens::Vector{Perm}
    one::Perm
    sims::Union{Nothing, NamedTuple}
    SimsGp(gens, one) = new(gens, one, nothing)
end

is_trivial(group::SimsGp) = all(isidentity, group.gens)

last_moved(group::SimsGp) = is_trivial(group) ? 0 : max(last_moved.(group.gens)...)

"""
    in(a, group)

Test whether permutation `a` belongs to `group` by sifting through the
Schreier-Sims stabilizer chain.
"""
function in(a::Perm, group::SimsGp)
    is_trivial(group) && return isidentity(a)
    s = sims(group)
    pos = findfirst(==(s.list[1]^a), s.list)
    isnothing(pos) && return false
    return a/s.reps[pos] in s.stab
end

"""
    closure(group, a)

Return the smallest `SimsGp` containing both `group` and the permutation `a`.
"""
function closure(group::SimsGp, a::Perm)
    a in group && return group
    return SimsGp(vcat(group.gens, a), group.one)
end

"""
    orbit_sims(aaa, x, under=^)

Compute the orbit of point `x` under generators `aaa` with action `under`,
recording a right transversal and the point stabilizer as a `SimsGp`.
Returns a named tuple with fields `list`, `reps`, and `stab`.
"""
function orbit_sims(aaa, x, under=^)
    list = [x]
    index = Dict(x => 1)
    reps = [aaa[1]^0]
    stab = SimsGp([], aaa[1]^0)
    for (i, y) in enumerate(list)
        for a in aaa
            z = under(y, a)
            t = onRight(reps[i], a)
            l = get!(index, z) do
                push!(list, z)
                push!(reps, t)
                length(list)
            end   # x^(reps[i] * a) = x^reps[l]
            t == reps[l] || (stab = closure(stab, t / reps[l]))
        end
    end
    return (list = list, reps = reps, stab = stab)
end

function sims(group::SimsGp)
    if isnothing(group.sims)
        group.sims = orbit_sims(group.gens, last_moved(group))
    end
    return group.sims
end

"""
    size(group)

Return the order of `group` as a `BigInt`, computed as the product of orbit
lengths along the stabilizer chain.
"""
function size(group::SimsGp)::BigInt
    is_trivial(group) && return 1
    length(sims(group).list) * size(sims(group).stab)
end

"""
    rand(group)

Return a uniformly random element of `group` by multiplying a random stabilizer
element with a random coset representative.
"""
function rand(group::SimsGp)
    is_trivial(group) && return group.one
    rand(sims(group).stab) * rand(sims(group).reps)
end

# Search the coset Gα·t for an element in H, using Gα's stabilizer chain.
function find_coset_rep_in(Gα::SimsGp, t::Perm, H::SimsGp)
    is_trivial(Gα) && return t ∈ H ? t : nothing
    s = sims(Gα)
    for v in s.reps
        c = find_coset_rep_in(s.stab, v * t, H)
        isnothing(c) || return c
    end
    return nothing
end

"""
    intersect(G::SimsGp, H::SimsGp)
    G ∩ H

Return generators of G ∩ H as a `SimsGp`, using an orbit-stabilizer recursion
with Schreier-Sims membership testing.

Uses `orbit_sims` to compute the orbit and stabilizer of a base point `α` in
both groups. For each `β` in `α^G ∩ α^H`, searches the coset `G_α · t_β` for
an element lying in H using `find_coset_rep_in`, which recurses through G_α's
own stabilizer chain rather than enumerating its elements. The stabilizer part
`G_α ∩ H_α` is found by the same algorithm recursively.

# Examples
```jldoctest
julia> using OrbitAl

julia> s = Perm([2,1,3,4]); t = Perm([1,3,2,4]); u = Perm([1,2,4,3]);

julia> G = SimsGp([s, t], one(s));  # S3 on {1,2,3}, order 6

julia> H = SimsGp([t, u], one(s));  # S3 on {2,3,4}, order 6

julia> K = G ∩ H;

julia> size(K)
2
```
"""
function intersect(G::SimsGp, H::SimsGp)
    is_trivial(G) && return SimsGp([], G.one)
    is_trivial(H) && return SimsGp([], G.one)

    α = max(last_moved(G), last_moved(H))

    oG = orbit_sims(G.gens, α)
    oH = orbit_sims(H.gens, α)

    H_orbit = Set(oH.list)

    coset_reps = Perm[]
    for i in eachindex(oG.list)
        β = oG.list[i]
        β == α && continue
        β ∈ H_orbit || continue
        c = find_coset_rep_in(oG.stab, oG.reps[i], H)
        isnothing(c) || push!(coset_reps, c)
    end

    K = intersect(oG.stab, oH.stab)

    gens = filter(!isidentity, vcat(K.gens, coset_reps))
    return SimsGp(gens, G.one)
end

## test data: Rubik's cube group
##
##                 +--------------+
##                 |              |
##                 |  1    2    3 |
##                 |              |
##                 |  4  top    5 |
##                 |              |
##                 |  6    7    8 |
##                 |              |
##  +--------------+--------------+--------------+--------------+
##  |              |              |              |              |
##  |  9   10   11 | 17   18   19 | 25   26   27 | 33   34   35 |
##  |              |              |              |              |
##  | 12  left  13 | 20 front  21 | 28 right  29 | 36  rear  37 |
##  |              |              |              |              |
##  | 14   15   16 | 22   23   24 | 30   31   32 | 38   39   40 |
##  |              |              |              |              |
##  +--------------+--------------+--------------+--------------+
##                 |              |
##                 | 41   42   43 |
##                 |              |
##                 | 44 bottom 45 |
##                 |              |
##                 | 46   47   48 |
##                 |              |
##                 +--------------+
##
cube = SimsGp([Perm(48, cycles) for cycles in [
  [[ 1, 3, 8, 6],[ 2, 5, 7, 4],[ 9,33,25,17],[10,34,26,18],[11,35,27,19]],
  [[ 9,11,16,14],[10,13,15,12],[ 1,17,41,40],[ 4,20,44,37],[ 6,22,46,35]],
  [[17,19,24,22],[18,21,23,20],[ 6,25,43,16],[ 7,28,42,13],[ 8,30,41,11]],
  [[25,27,32,30],[26,29,31,28],[ 3,38,43,19],[ 5,36,45,21],[ 8,33,48,24]],
  [[33,35,40,38],[34,37,39,36],[ 3, 9,46,32],[ 2,12,47,29],[ 1,14,48,27]],
  [[41,43,48,46],[42,45,47,44],[14,22,30,38],[15,23,31,39],[16,24,32,40]],
]], Perm(48))

end # module
