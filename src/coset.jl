#############################################################################
##
#A  coset.jl                                                    orbits-julia
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  Action on the cosets of a subgroup
##
module coset

import Base: ==, *
import .Main.PermGp
import permutation: Perm

export Coset
export cosets

struct Coset
    subgp::PermGp
    element::Perm
end

*(coset::Coset, perm::Perm) = Coset(coset.subgp, coset.element * perm)

==(coset::Coset, other::Coset) =
    coset.subgp == other.subgp && coset.element/other.element in coset.subgp

cosets(group::PermGp, subgp::PermGp) =
    orbit(group.gens, Coset(subgp, group.one), onRight)

end # module
