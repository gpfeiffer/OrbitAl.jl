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

abstract type APermGp end

## Perm group
struct PermGp <: APermGp
    gens::Vector{Perm}
    one::Perm
end

## elements
elements(group::APermGp) = sort(orbit(group.gens, group.one, onRight))

## equality, comparison
==(group::APermGp, other::APermGp) = elements(group) == elements(other)
isless(group::APermGp, other::APermGp) = elements(group) < elements(other)

##  closure
closure(group::APermGp, a::Perm) = PermGp(union(group.gens, [a]), group.one)
onGroups(x, a) = closure(x, a)

##  subgroups
subgroups(group) = orbit(elements(group), PermGp([], group.one), onGroups)

## conjugation
class(group::APermGp, a::Perm) = orbit(group.gens, a, onPoints)
^(a::Perm, group::APermGp) = Orbit(group, sort(class(group, a)))

^(group::APermGp, a::Perm) = PermGp([x^a for x in group.gens], group.one)

subgpClass(gp::APermGp, subgp::APermGp) = orbit(gp.gens, subgp, onPoints)
^(subgp::APermGp, group::APermGp) = Orbit(group, sort(subgpClass(group, subgp)))

## conjugacy classes
onClasses(x, a) = (x.elts[1] * a)^(x.group)

function conjClasses(gp::APermGp)
    orbit(orbitx(gp.gens, copy(gp.gens), onPoints), gp.one^gp, onClasses)
end

onSubgpClasses(x, a) = onGroups(x.elts[1], a)^(x.group)

function subgpClasses(gp::APermGp)
    orbit(elements(gp), PermGp([], gp.one)^gp, onSubgpClasses)
end

## is trivial?
is_trivial(group::APermGp) = all(isidentity, group.gens)

## largest moved point
last_moved(group::APermGp) = max(last_moved.(group.gens)...)

## size of group
function sizeOfGroup(group::APermGp)
    is_trivial(group) && return 1
    x = last_moved(group)
    orb = orbit_with_stabilizer(group.gens, x, onPoints)
    stab = PermGp(setdiff(orb.stab, [group.one]), group.one)
    return sizeOfGroup(stab) * length(orb.reps)
end
size(group::APermGp) = sizeOfGroup(group)

## random group element
function randomGroupElement(group::APermGp)
    is_trivial(group) && return group.one
    x = last_moved(group)
    orb = orbit_with_stabilizer(group.gens, x, onPoints)
    stab = PermGp(setdiff(orb.stab, [group.one]), group.one)
    return randomGroupElement(stab) * rand(orb.reps)
end
rand(group::APermGp) = randomGroupElement(group)

## membership test
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

end # module
