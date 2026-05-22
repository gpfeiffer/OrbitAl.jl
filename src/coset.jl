#############################################################################
##
#A  coset.jl                                                         OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  Action on the cosets of a subgroup
##
module coset

import Base: ==, *

using ..permutation
using ..orbits
using ..permgroup

import ..permutation: Perm
import ..permgroup: APermGp

export Coset, cosets

"""
    Coset(subgp, element)

A right coset of `subgp` represented by `element`.
Two cosets are equal if their representatives differ by an element of `subgp`.
"""
struct Coset
    subgp::APermGp
    element::Perm
end

*(c::Coset, a::Perm) = Coset(c.subgp, c.element * a)

==(c::Coset, other::Coset) =
    c.subgp == other.subgp && c.element / other.element ∈ c.subgp

"""
    cosets(group, subgp)

Return the right cosets of `subgp` in `group` as a vector of `Coset` objects,
computed by orbit enumeration under right multiplication by `group`'s generators.
"""
cosets(group::APermGp, subgp::APermGp) =
    orbit(group.gens, Coset(subgp, group.one), onRight)

end # module
