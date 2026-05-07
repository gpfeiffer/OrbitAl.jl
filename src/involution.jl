#############################################################################
##
#A  involution.jl                                                     OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  Action on involutions and classes of involutions
##
module involution

using ..orbits
import ..coxeter: reflections

export onInvolutions, involutions, onInvolutionClasses, involutionClasses

"""
    onInvolutions(x, s)

The involution action of generator `s` on a group element `x`:
- If `x` and `s` commute (`x*s == s*x`): return `x*s` (right multiplication).
- Otherwise: return `s⁻¹*x*s` (conjugation by `s`).

Both cases preserve involutions: if `x² = 1` then the result squares to 1
in either case. The orbit of the identity under this action via `involutions`
yields all involutions of the group.

# Examples
```jldoctest
julia> using OrbitAl

julia> W = CoxeterGp(cartanMat("A", 2));

julia> s = W.gens[1];

julia> onInvolutions(W.one, s) == s
true
```
"""
onInvolutions(x, s) = x * s == s * x ? onRight(x, s) : x^s

"""
    involutions(W)

Return all involutions in `W` (elements `x` with `x² = 1`) as the orbit of
the identity under the involution action `onInvolutions`.

# Examples
```jldoctest
julia> using OrbitAl

julia> W = CoxeterGp(cartanMat("A", 2));

julia> length(involutions(W))
4
```
"""
involutions(W) = orbit(W.gens, W.one, onInvolutions)

"""
    onInvolutionClasses(x, a)

The action of reflection `a` on a conjugacy class of involutions `x`
(an `Orbit`). Let `y` be the canonical (minimum) representative of `x`:
- If `a` centralises `y` (`y^a == y`, equivalently `y*a == a*y`): `y*a` is a
  new involution, return its conjugacy class.
- Otherwise: `a` conjugates `y` within its class, so the class is fixed;
  return `x` unchanged.

The orbit of the identity class under this action via `involutionClasses`
yields all conjugacy classes of involutions.

# Examples
```jldoctest
julia> using OrbitAl

julia> W = CoxeterGp(cartanMat("A", 2));

julia> cls = W.one^W;

julia> size(onInvolutionClasses(cls, reflections(W)[1]))
3
```
"""
function onInvolutionClasses(x, a)
    y = x.elts[1]
    # y^a == y means a centralises y, so y*a is a new involution
    y^a != y ? x : onRight(y, a)^x.group
end

"""
    involutionClasses(W)

Return all conjugacy classes of involutions in `W` as the orbit of the
identity class under the action `onInvolutionClasses` by reflections.
Each element of the result is an `Orbit` (a conjugacy class).

# Examples
```jldoctest
julia> using OrbitAl

julia> W = CoxeterGp(cartanMat("A", 2));

julia> length(involutionClasses(W))
2
```
"""
involutionClasses(W) = orbit(reflections(W), W.one^W, onInvolutionClasses)

end # module
