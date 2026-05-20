#############################################################################
##
#A  permgroup.jl                                                      OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  Lightweight Permutation Groups with Applications of Orbit Algorithms
##
module permgroup

using ..permutation
using ..orbits

import Base: in, isless, rand, size, ==, ^
import ..permutation: isidentity, last_moved

export APermGp, PermGp
export elements, conjClasses, closure, subgroups, subgpClasses
export sizeOfGroup, randomGroupElement, memberOfGroup
export isPrimePower, zuppos

"""
    APermGp

Abstract supertype for all permutation group types.
"""
abstract type APermGp end

"""
    PermGp(gens, one)

A permutation group with generators `gens::Vector{Perm}` and identity element `one`.

# Examples
```jldoctest
julia> using OrbitAl

julia> s = Perm([2,1,3]); t = Perm([1,3,2]);

julia> G = PermGp([s, t], one(s));

julia> sizeOfGroup(G)
6
```
"""
struct PermGp <: APermGp
    gens::Vector{Perm}
    one::Perm
end

"""
    elements(group)

Return all elements of `group` as a sorted vector.
"""
elements(group::APermGp) = sort(orbit(group.gens, group.one, onRight))

## equality, comparison
==(group::APermGp, other::APermGp) = elements(group) == elements(other)
isless(group::APermGp, other::APermGp) = elements(group) < elements(other)

"""
    closure(group, a)

Return a new `PermGp` formed by adding generator `a` to `group`.
"""
closure(group::APermGp, a::Perm) = PermGp(union(group.gens, [a]), group.one)
onGroups(x, a) = closure(x, a)

"""
    subgroups(group)

Return all subgroups of `group` as a vector of `PermGp` values,
computed by orbit enumeration under the closure action.
"""
subgroups(group) = orbit(elements(group), PermGp([], group.one), onGroups)

## conjugation
class(group::APermGp, a::Perm) = orbit(group.gens, a, onPoints)
^(a::Perm, group::APermGp) = Orbit(group, sort(class(group, a)))

^(group::APermGp, a::Perm) = PermGp([x^a for x in group.gens], group.one)

subgpClass(gp::APermGp, subgp::APermGp) = orbit(gp.gens, subgp, onPoints)
^(subgp::APermGp, group::APermGp) = Orbit(group, sort(subgpClass(group, subgp)))

## conjugacy classes
onClasses(x, a) = (x.elts[1] * a)^(x.group)

"""
    conjClasses(gp)

Return the conjugacy classes of `gp` as a vector of `Orbit` objects,
each containing the elements of one class.
"""
function conjClasses(gp::APermGp)
    orbit(orbitx(gp.gens, copy(gp.gens), onPoints), gp.one^gp, onClasses)
end

onSubgpClasses(x, a) = onGroups(x.elts[1], a)^(x.group)

"""
    subgpClasses(gp)

Return the conjugacy classes of subgroups of `gp` as a vector of `Orbit` objects.
"""
function subgpClasses(gp::APermGp)
    orbit(elements(gp), PermGp([], gp.one)^gp, onSubgpClasses)
end

## is trivial?
is_trivial(group::APermGp) = all(isidentity, group.gens)

## largest moved point
last_moved(group::APermGp) = max(last_moved.(group.gens)...)

"""
    sizeOfGroup(group)

Return the order of `group` using the orbit-stabilizer theorem,
without enumerating all elements.

# Examples
```jldoctest
julia> using OrbitAl

julia> s = Perm([2,1,3]); t = Perm([1,3,2]);

julia> sizeOfGroup(PermGp([s, t], one(s)))
6
```
"""
function sizeOfGroup(group::APermGp)
    is_trivial(group) && return 1
    x = last_moved(group)
    orb = orbit_with_stabilizer(group.gens, x, onPoints)
    stab = PermGp(setdiff(orb.stab, [group.one]), group.one)
    return sizeOfGroup(stab) * length(orb.reps)
end
size(group::APermGp) = sizeOfGroup(group)

"""
    randomGroupElement(group)

Return a uniformly random element of `group` using a recursive
orbit-stabilizer decomposition.
"""
function randomGroupElement(group::APermGp)
    is_trivial(group) && return group.one
    x = last_moved(group)
    orb = orbit_with_stabilizer(group.gens, x, onPoints)
    stab = PermGp(setdiff(orb.stab, [group.one]), group.one)
    return randomGroupElement(stab) * rand(orb.reps)
end
rand(group::APermGp) = randomGroupElement(group)

"""
    memberOfGroup(group, perm)

Return `true` if `perm` belongs to `group`, using a recursive sifting algorithm.

# Examples
```jldoctest
julia> using OrbitAl

julia> s = Perm([2,1,3]); t = Perm([1,3,2]);

julia> G = PermGp([s, t], one(s));

julia> memberOfGroup(G, s)
true

julia> memberOfGroup(G, Perm([2,3,4,1]))
false
```
"""
function memberOfGroup(group::APermGp, perm::Perm)
    is_trivial(group) && return perm == group.one
    x = last_moved(group)
    orb = orbit_with_stabilizer(group.gens, x, onPoints)
    l = findfirst(==(x^perm), orb.list)
    isnothing(l) && return false
    stab = PermGp(setdiff(orb.stab, [group.one]), group.one)
    return memberOfGroup(stab, perm / orb.reps[l])
end
in(a::Perm, group::APermGp) = memberOfGroup(group, a)

"""
    isPrimePower(n)

Return `true` if `n > 1` is a power of a single prime, i.e., `n = p^k` for some prime `p`
and `k ≥ 1`.
"""
function isPrimePower(n::Int)
    n > 1 || return false
    p = first(p for p in 2:n if n % p == 0)
    while n % p == 0; n ÷= p; end
    return n == 1
end

"""
    zuppos(group)

Return one generator for each cyclic subgroup of prime-power order (zuppo) of `group`.
The returned generators suffice to find all subgroups of `group` under the closure action.
"""
function zuppos(group::APermGp)
    list = [group.one]
    Z = Perm[]
    for g in list
        for a in group.gens
            z = g * a
            z ∈ list || push!(list, z)
        end
        g == group.one && continue
        ord = order(g)
        isPrimePower(ord) || continue
        any(order(z) == ord && g ∈ orbit([z], group.one, onRight) for z in Z) && continue
        push!(Z, g)
    end
    return Z
end
end # module
