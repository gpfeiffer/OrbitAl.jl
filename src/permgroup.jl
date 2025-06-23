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

export PermGp, elements, conjClasses, closure, subgroups, subgpClasses
export sizeOfGroup, randomGroupElement, memberOfGroup

## Perm group
struct PermGp
    gens::Vector{Perm}
    one::Perm
end

## elements
elements(group::PermGp) = sort(orbit(group.gens, group.one, onRight))

## equality, comparison
==(group::PermGp, other::PermGp) = elements(group) == elements(other)
isless(group::PermGp, other::PermGp) = elements(group) < elements(other)

##  closure
closure(group::PermGp, a::Perm) = PermGp(union(group.gens, [a]), group.one)
onGroups(x, a) = closure(x, a)

##  subgroups
subgroups(group) = orbit(elements(group), PermGp([], group.one), onGroups)

## conjugation
class(group::PermGp, a::Perm) = orbit(group.gens, a, onPoints)
^(a::Perm, group::PermGp) = Orbit(group, sort(class(group, a)))

^(group::PermGp, a::Perm) = PermGp([x^a for x in group.gens], group.one)

subgpClass(gp::PermGp, subgp::PermGp) = orbit(gp.gens, subgp, onPoints)
^(subgp::PermGp, group::PermGp) = Orbit(group, sort(subgpClass(group, subgp)))

## conjugacy classes
onClasses(x, a) = (x.elts[1] * a)^(x.group)

function conjClasses(gp::PermGp)
    orbit(orbitx(gp.gens, gp.gens, onPoints), gp.one^gp, onClasses)
end

onSubgpClasses(x, a) = onGroups(x.elts[1], a)^(x.group)

function subgpClasses(gp::PermGp)
    orbit(elements(gp), PermGp([], gp.one)^gp, onSubgpClasses)
end

## is trivial?
is_trivial(group::PermGp) = all(isidentity, group.gens)

## largest moved point
last_moved(group::PermGp) = max(last_moved.(group.gens)...)

## size of group
function sizeOfGroup(group::PermGp)
    is_trivial(group) && return 1
    x = last_moved(group)
    orb = orbit_with_stabilizer(group.gens, x, onPoints)
    stab = PermGp(setdiff(orb.stab, [group.one]), group.one)
    return sizeOfGroup(stab) * length(orb.reps)
end
size(group::PermGp) = sizeOfGroup(group)

## random group element
function randomGroupElement(group::PermGp)
    is_trivial(group) && return group.one
    x = last_moved(group)
    orb = orbit_with_stabilizer(group.gens, x, onPoints)
    stab = PermGp(setdiff(orb.stab, [group.one]), group.one)
    return randomGroupElement(stab) * rand(orb.reps)
end
rand(group::PermGp) = randomGroupElement(group)

## membership test
function memberOfGroup(group::PermGp, perm::Perm)
    is_trivial(group) && return perm == group.one
    x = last_moved(group)
    orb = orbit_with_stabilizer(group.gens, x, onPoints)
    l = findfirst(==(x^perm), orb.list)
    isnothing(l) && return false
    stab = PermGp(setdiff(orb.stab, [group.one]), group.one)
    return memberOfGroup(stab, perm / orb.reps[l])
end
in(a::Perm, group::PermGp) = memberOfGroup(a, group)

end # module
