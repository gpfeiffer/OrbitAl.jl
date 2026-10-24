#############################################################################
##
#A  palg.jl                                                           OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  Algebras on a prefix-closed basis of words
##
##  An algebra `A`, free of rank `d` over a commutative ring `Lambda`, given
##  by generators `p_1, .., p_r` and a basis of positive words in them,
##
##      b_1 = 1,  ..,  b_d,   b_x = b_(parent[x]) T_(letter[x]),
##
##  closed under prefixes.  What is stored is the sparse multiplication by one
##  generator on each side, `right[j][x] = b_x T_j` and `left[j][x] = T_j b_x`;
##  every other product is built from those.  The algebra need not commute,
##  so the two tables are different data.
##
##  Each generator has a minimal polynomial, `params[j] = [c_0, .., c_(k-1)]`
##  with the constant term first:
##
##      T_j^k = c_(k-1) T_j^(k-1) + .. + c_1 T_j + c_0,
##
##  and `c_0` is a unit, which is what makes `T_j^-1` an element of `A`.
##  Generators that are to share a minimal polynomial share a vector.
##
##  `cyclicAlg(params)` is the one-generator case, where the basis is
##  `1, u, .., u^(d-1)` and the algebra is commutative.  `ugen`, `uinv`,
##  `upow` and `squareParams` belong to it alone and assert `iscyclic`.
##
##  The coefficient type `T` carries the ground ring: `Zn{p}` for a run over
##  `F_p`, a polynomial type for a symbolic one, and nothing here branches on
##  which.  Zero and one are taken from a parameter instance, not from
##  `zero(T)`, so a type that defines them only on instances still works.
##
##  This is the Hecke algebra of a parabolic subgroup when the generators are
##  braid generators, but nothing below assumes that.  Unlike
##  [`OrbitAl.hecke`](@ref), which stores no elements and so serves a Coxeter
##  group as large as `E_8`, a `PAlg` stores both multiplication tables and is
##  for an algebra whose basis can be held: ranks into the tens of thousands.
##
module palg

export PAlg, PElt
export basis, element, basisInverse, tau, coefficient, mul, scale, inverse
export coeffs, fromcoeffs, multMatrix, zeroc, onec
export cyclicAlg, iscyclic, ugen, uinv, upow, squareParams
export withAcc, addall, sumall
export rmulgen, lmulgen, mulinv, rmulword, lmulword
export findUnitInverse

"""
    PElt{T}

An element of a [`PAlg`](@ref), `sum val[k] b_(idx[k])` with `idx` strictly
increasing and `val` free of zeros.  That form is canonical, so equality is
structural.
"""
struct PElt{T}
    idx::Vector{Int32}
    val::Vector{T}
end

Base.:(==)(x::PElt, y::PElt) = x.idx == y.idx && x.val == y.val
Base.hash(x::PElt, h::UInt) = hash(x.val, hash(x.idx, h))
Base.iszero(x::PElt) = isempty(x.idx)

function Base.show(io::IO, x::PElt)
    isempty(x.idx) && return print(io, "0")
    cf(c) = isone(c) ? "" : "(" * sprint(show, c; context = :typeinfo => typeof(c)) * ")"
    print(io, join([cf(c) * "b$i" for (i, c) in zip(x.idx, x.val)], " + "))
end

"""
    coefficient(x, i, z)

The coefficient of `b_i` in `x`, which is total: zero off the support.  The
zero is the one `z` supplies, since `T` need not define `zero(::Type{T})`.
"""
function coefficient(x::PElt{T}, i::Integer, z::T) where T
    k = searchsortedfirst(x.idx, Int32(i))
    (k <= length(x.idx) && x.idx[k] == i) ? x.val[k] : z
end

# c_0^-1, where an integer coefficient is a unit only at +-1
inverse(x::Integer) = abs(x) == 1 ? x : throw(ArgumentError("$x is not a unit"))
inverse(x) = inv(x)

##  c x for a scalar c of Lambda, which is central, so there is no side to
##  choose.  Lambda can have zero divisors, so a scaled entry may vanish.
function scale(c, x::PElt{T}) where T
    isone(c) && return x
    idx = Int32[]; val = T[]
    for k in eachindex(x.idx)
        v = c * x.val[k]
        iszero(v) || (push!(idx, x.idx[k]); push!(val, v))
    end
    PElt{T}(idx, val)
end

##  x + sgn y, sgn = ±1, merging the two supports in one pass.  A sum that
##  cancels is dropped, so the result is canonical.
function lincomb(x::PElt{T}, y::PElt{T}, sgn::Int) where T
    idx = Int32[]; val = T[]
    i = 1; j = 1
    nx = length(x.idx); ny = length(y.idx)
    @inbounds while i <= nx || j <= ny
        if j > ny || (i <= nx && x.idx[i] < y.idx[j])
            push!(idx, x.idx[i]); push!(val, x.val[i])
            i += 1
        elseif i > nx || y.idx[j] < x.idx[i]
            push!(idx, y.idx[j]); push!(val, sgn == 1 ? y.val[j] : -y.val[j])
            j += 1
        else
            v = sgn == 1 ? x.val[i] + y.val[j] : x.val[i] - y.val[j]
            iszero(v) || (push!(idx, x.idx[i]); push!(val, v))
            i += 1; j += 1
        end
    end
    PElt{T}(idx, val)
end

Base.:+(x::PElt{T}, y::PElt{T}) where T =
    iszero(y) ? x : iszero(x) ? y : lincomb(x, y, 1)
Base.:-(x::PElt{T}, y::PElt{T}) where T = iszero(y) ? x : lincomb(x, y, -1)
Base.:-(x::PElt{T}) where T = PElt{T}(x.idx, T[-v for v in x.val])
Base.:*(c::Number, x::PElt) = scale(c, x)
Base.:*(x::PElt, c::Number) = scale(c, x)
Base.:*(c::T, x::PElt{T}) where T = scale(c, x)
Base.:*(x::PElt{T}, c::T) where T = scale(c, x)

#############################################################################
##
##  The accumulator.
##
##  A product reaches each basis element many times, so the terms are summed
##  on a dense buffer and turned into an element once, at the end.  Only the
##  slots touched are visited, so the cost is the support and not the rank.
##
##  The buffers are POOLED and the pool is indexed by depth, because a product
##  may be computed while a caller holds one: `mul` multiplies along a memo of
##  prefixes, and a caller may nest a walk inside a walk.  One shared buffer
##  would be silently corrupted by that, and the result would look plausible.
##  Nothing here is threaded and the pool does not pretend to be thread-safe.
##
mutable struct Acc{T}
    buf::Vector{T}
    touched::Vector{Int32}
    mark::BitVector
end

Acc{T}(d::Integer, z::T) where T = Acc{T}(fill(z, d), Int32[], falses(d))

##  add c s into slot y, recording the slot the first time it is reached
@inline function addmul!(a::Acc{T}, y::Integer, c::T, s::T) where T
    a.mark[y] || (a.mark[y] = true; push!(a.touched, y))
    a.buf[y] += c * s
    nothing
end

##  empty the touched slots into a canonical element, leaving the buffer zero
function collect!(a::Acc{T}, z::T) where T
    sort!(a.touched)
    idx = Int32[]; val = T[]
    for y in a.touched
        v = a.buf[y]
        a.buf[y] = z
        a.mark[y] = false
        iszero(v) || (push!(idx, y); push!(val, v))
    end
    empty!(a.touched)
    PElt{T}(idx, val)
end

#############################################################################

"""
    PAlg{T}

An algebra on a prefix-closed basis of positive words, with the right and
left multiplications by each generator stored.

`params[j] = [c_0, .., c_(k-1)]`, constant term first, is the minimal
polynomial of the generator `p_j`: `T^k = c_(k-1) T^(k-1) + .. + c_1 T + c_0`.
Generators that share a minimal polynomial share a vector.
"""
mutable struct PAlg{T}
    name::String
    d::Int                                  # the rank over Lambda
    ngens::Int                              # the generators
    params::Vector{Vector{T}}               # one minimal polynomial per generator
    words::Vector{Vector{Int}}              # b_x as a positive word in 1 .. ngens
    parent::Vector{Int}
    letter::Vector{Int}
    right::Vector{Vector{PElt{T}}}          # right[j][x] = b_x T_j
    left::Vector{Vector{PElt{T}}}           # left[j][x]  = T_j b_x
    suf::Vector{Int}                        # b_x's word as a node of the trie of all suffixes
    sufnext::Vector{Int}                    # the node of `s w` maps to the node of `w`; node 1 is the empty word
    sufletter::Vector{Int}
    bylength::Vector{Int}                   # the basis, by increasing word length
    accs::Vector{Acc{T}}                    # the buffer pool, indexed by depth
    depth::Int
    cyclic::Bool
end

Base.show(io::IO, A::PAlg) = print(io, A.name)

##  the zero and one of Lambda, from a parameter instance: a coefficient type
##  that defines them only on instances still works
zeroc(A::PAlg) = zero(A.params[1][1])
onec(A::PAlg) = one(A.params[1][1])

Base.zero(A::PAlg{T}) where T = PElt{T}(Int32[], T[])
Base.one(A::PAlg{T}) where T = PElt{T}(Int32[1], T[onec(A)])

"""
    basis(A, x)

The basis element `b_x`.
"""
basis(A::PAlg{T}, x::Integer) where T = PElt{T}(Int32[x], T[onec(A)])

"""
    tau(A, h)

`tau(h)`, the coefficient of `b_1 = 1`: the canonical trace of `A`.
"""
tau(A::PAlg{T}, h::PElt{T}) where T = coefficient(h, 1, zeroc(A))

"""
    coefficient(A, h, i)

The coefficient of `b_i` in `h`, with the zero taken from `A`.
"""
coefficient(A::PAlg{T}, h::PElt{T}, i::Integer) where T = coefficient(h, i, zeroc(A))

"""
    iscyclic(A)

Whether `A` has one generator, so that it is commutative and `ugen`, `uinv`,
`upow` and `squareParams` apply.
"""
iscyclic(A::PAlg) = A.cyclic

"""
    PAlg(name, params, words, right, left)

The algebra on `words`, which must be closed under prefixes and ordered so
that a word comes after its prefix.
"""
function PAlg(name::AbstractString, params::Vector{Vector{T}},
              words::Vector{Vector{Int}}, right, left) where T
    d = length(words)
    d >= 1 || throw(ArgumentError("the basis is empty"))
    isempty(words[1]) || throw(ArgumentError("the first basis element is not 1"))
    index = Dict(w => x for (x, w) in enumerate(words))
    length(index) == d || throw(ArgumentError("the basis repeats a word"))
    parent = [isempty(w) ? 0 : get(index, w[1:end-1]) do
                  throw(ArgumentError("the basis is not closed under prefixes: $w"))
              end for w in words]
    all(x -> x == 1 || parent[x] < x, 1:d) ||
        throw(ArgumentError("the basis is not ordered along its tree"))

    ##  the trie of all suffixes, so that a product can be expanded from
    ##  whichever side is cheaper
    node = Dict(Int[] => 1)
    sufnext = [0]
    sufletter = [0]
    function sufnode(w)
        get!(node, w) do
            t = sufnode(w[2:end])
            push!(sufnext, t)
            push!(sufletter, w[1])
            length(sufnext)
        end
    end
    suf = [sufnode(w) for w in words]

    PAlg{T}(name, d, length(right), params, words, parent,
            [isempty(w) ? 0 : w[end] for w in words], right, left,
            suf, sufnext, sufletter, sortperm(length.(words)),
            [Acc{T}(d, zero(params[1][1]))], 0,
            length(right) == 1 && d == length(params[1]))
end

"""
    cyclicAlg(params)

The algebra of one generator `u`: the basis is `1, u, .., u^(d-1)` and the
algebra is commutative, so the two multiplication tables coincide.
`params = [c_0, .., c_(d-1)]`, constant term first, and `c_0` must be a unit
since `u^-1` is read off it.
"""
function cyclicAlg(params::AbstractVector{T}) where T
    k = length(params)
    k >= 1 || throw(ArgumentError("params is empty: deg must be >= 1"))
    iszero(params[1]) && throw(ArgumentError("c_0 = 0: u is not a unit"))
    o = one(params[1])
    words = [fill(1, j) for j in 0:k-1]
    ##  b_(j+1) T_1 = u^(j+1), which is a basis element until it folds back
    rows = [j < k - 1 ? PElt{T}(Int32[j + 2], T[o]) : PElt{T}(Int32.(1:k), collect(T, params))
            for j in 0:k-1]
    PAlg("Z$k", [collect(T, params)], words, [rows], [rows])
end

#############################################################################
##
##  The product.
##
"""
    withAcc(f, A)

Run `f` on a buffer of its own, so that a product computed inside another one
does not corrupt it.
"""
function withAcc(f, A::PAlg{T}) where T
    A.depth += 1
    A.depth <= length(A.accs) || push!(A.accs, Acc{T}(A.d, zeroc(A)))
    try
        f(A.accs[A.depth])
    finally
        A.depth -= 1
    end
end

##  x T_j with mats = A.right[j], or T_j x with mats = A.left[j]
function mulgen(A::PAlg{T}, x::PElt{T}, mats::Vector{PElt{T}}) where T
    withAcc(A) do a
        for k in eachindex(x.idx)
            c = x.val[k]
            row = mats[x.idx[k]]
            for j in eachindex(row.idx)
                addmul!(a, row.idx[j], c, row.val[j])
            end
        end
        collect!(a, zeroc(A))
    end
end

"""
    rmulgen(A, x, j)

`x T_j`.
"""
rmulgen(A::PAlg, x::PElt, j::Int) = mulgen(A, x, A.right[j])

"""
    lmulgen(A, j, x)

`T_j x`.
"""
lmulgen(A::PAlg, j::Int, x::PElt) = mulgen(A, x, A.left[j])

"""
    mulinv(A, x, mats, j)

`x T_j^-1` with `mats = A.right[j]`, or `T_j^-1 x` with `mats = A.left[j]`,
by Horner on

    c_0 T^-1 = T^(k-1) - c_(k-1) T^(k-2) - ... - c_1.

`c_0` is the only coefficient that has to be a unit.
"""
function mulinv(A::PAlg{T}, x::PElt{T}, mats::Vector{PElt{T}}, j::Int) where T
    c = A.params[j]
    img = x
    for i in length(c)-1:-1:1
        img = mulgen(A, img, mats) - scale(c[i+1], x)
    end
    scale(inverse(c[1]), img)
end

##  A letter is a signed label: -s is the inverse of the generator s.  A word
##  multiplies on the right letter by letter from the left, and on the left
##  letter by letter from the right, so both read in the order written.
rmulletter(A::PAlg, x::PElt, s::Int) = s > 0 ? rmulgen(A, x, s) : mulinv(A, x, A.right[-s], -s)
lmulletter(A::PAlg, s::Int, x::PElt) = s > 0 ? lmulgen(A, s, x) : mulinv(A, x, A.left[-s], -s)

"""
    rmulword(A, x, w)

`x T_w` for a word `w` in the signed labels `±1 .. ±ngens`.
"""
rmulword(A::PAlg, x::PElt, w) = foldl((h, s) -> rmulletter(A, h, s), w; init = x)

"""
    lmulword(A, w, x)

`T_w x` for a word `w` in the signed labels `±1 .. ±ngens`.
"""
lmulword(A::PAlg, w, x::PElt) = foldr((s, h) -> lmulletter(A, s, h), w; init = x)

"""
    element(A, w)

`T_w` for a word `w` in the signed labels `±1 .. ±ngens`.
"""
element(A::PAlg, w) = rmulword(A, one(A), w)

"""
    basisInverse(A, x)

`b_x^-1 = T_(s_l)^-1 .. T_(s_1)^-1` for `b_x = T_(s_1) .. T_(s_l)`.  A basis
element is a product of generators, hence a unit.
"""
basisInverse(A::PAlg, x::Integer) =
    lmulword(A, [-s for s in Iterators.reverse(A.words[x])], one(A))

"""
    addall(A, terms)

`sum c_k x_k` for `terms` of scalar-element pairs `(c_k, x_k)`, on one pass
of the buffer rather than pairwise merges.
"""
function addall(A::PAlg{T}, terms) where T
    withAcc(A) do a
        for (c, x) in terms, k in eachindex(x.idx)
            addmul!(a, x.idx[k], c, x.val[k])
        end
        collect!(a, zeroc(A))
    end
end

"""
    sumall(A, xs)

The sum of a list of elements, on one pass of the buffer.
"""
function sumall(A::PAlg{T}, xs) where T
    withAcc(A) do a
        for x in xs, k in eachindex(x.idx)
            a.mark[x.idx[k]] || (a.mark[x.idx[k]] = true; push!(a.touched, x.idx[k]))
            a.buf[x.idx[k]] += x.val[k]
        end
        collect!(a, zeroc(A))
    end
end

##  A cyclic algebra is Lambda[u]/(minimal polynomial), where the product is
##  a convolution of the coordinate sequences followed by a fold through the
##  minimal polynomial, HIGHEST POWER FIRST: folding u^i adds into slots
##  i-d .. i-1, all of them lower, so one descending pass reduces everything.
##  That is cheaper than the general memo at the ranks the cyclic case uses.
function cyclicmul(A::PAlg{T}, x::PElt{T}, y::PElt{T}) where T
    d = A.d
    z = zeroc(A)
    c = A.params[1]
    p = fill(z, 2d - 1)
    for i in eachindex(x.idx), j in eachindex(y.idx)
        p[x.idx[i]+y.idx[j]-1] += x.val[i] * y.val[j]
    end
    for i in (2d-1):-1:(d+1)
        t = p[i]
        iszero(t) && continue
        for m in 1:d
            p[i-d+m-1] += c[m] * t
        end
    end
    fromcoeffs(A, p[1:d])
end

"""
    mul(A, x, y)

`x y`.  Expanded along the prefixes of the words of the support of `y`, or
along the suffixes of the words of the support of `x`, whichever is cheaper:
the memo then costs one generator multiplication per word reached rather than
one per pair of basis elements.
"""
function mul(A::PAlg{T}, x::PElt{T}, y::PElt{T}) where T
    (iszero(x) || iszero(y)) && return zero(A)
    length(y.idx) == 1 && y.idx[1] == 1 && return scale(y.val[1], x)
    length(x.idx) == 1 && x.idx[1] == 1 && return scale(x.val[1], y)
    A.cyclic && return cyclicmul(A, x, y)

    need = Set{Int32}()   # x b_w for the words of supp(y) and their prefixes
    for w in y.idx
        while w != 1 && !(w in need)
            push!(need, w)
            w = A.parent[w]
        end
    end
    sufs = Set{Int}()    # T_u y for the suffixes of the words of supp(x)
    for v in x.idx
        t = A.suf[v]
        while t != 1 && !(t in sufs)
            push!(sufs, t)
            t = A.sufnext[t]
        end
    end

    costR = length(need) * min(A.d, 8 * length(x.idx))
    costL = length(sufs) * min(A.d, 8 * length(y.idx))
    if costR <= costL
        memo = Dict{Int32,PElt{T}}(Int32(1) => x)
        for w in sort!(collect(need))
            memo[w] = rmulgen(A, memo[A.parent[w]], A.letter[w])
        end
        addall(A, [(y.val[k], memo[y.idx[k]]) for k in eachindex(y.idx)])
    else
        memo = Dict{Int,PElt{T}}(1 => y)
        left(t) = get!(() -> lmulgen(A, A.sufletter[t], left(A.sufnext[t])), memo, t)
        addall(A, [(x.val[k], left(A.suf[x.idx[k]])) for k in eachindex(x.idx)])
    end
end

#############################################################################
##
##  Coordinates.
##
"""
    coeffs(A, x)

The `d` coordinates of `x` on the basis, zeros included: the dense form, for
a caller that wants a flat coefficient array.
"""
function coeffs(A::PAlg{T}, x::PElt{T}) where T
    out = fill(zeroc(A), A.d)
    for k in eachindex(x.idx)
        out[x.idx[k]] = x.val[k]
    end
    out
end

"""
    fromcoeffs(A, l)

The element with coordinates `l`, which must have length `d`.
"""
function fromcoeffs(A::PAlg{T}, l::AbstractVector) where T
    length(l) == A.d ||
        throw(ArgumentError("$(length(l)) coordinates, expected $(A.d)"))
    idx = Int32[]
    val = T[]
    for k in 1:A.d
        iszero(l[k]) || (push!(idx, Int32(k)); push!(val, T(l[k])))
    end
    PElt{T}(idx, val)
end

"""
    multMatrix(A, x)

The matrix of multiplication by `x` on the basis: row `i` is `b_i x`.  Its
entries are in `Lambda`, which is commutative whatever `A` is, so the
determinant of this matrix is the norm of `x`.
"""
function multMatrix(A::PAlg{T}, x::PElt{T}) where T
    mat = Matrix{T}(undef, A.d, A.d)
    for i in 1:A.d
        mat[i, :] = coeffs(A, mul(A, basis(A, i), x))
    end
    mat
end

#############################################################################
##
##  Units.
##
##  A basis element is a unit, and so is a product of a basis element with the
##  inverse of another.  `findUnitInverse` recognises `c b`, `c b^-1 b'` and
##  `c b b'^-1` for a unit `c` of `Lambda`, and verifies every candidate by
##  multiplying it out; it is the test a search for an invertible pivot wants,
##  where the norm through `multMatrix` is correct and too expensive.
##
##  Which coefficients are units of `Lambda` depends on the ground ring, so
##  `isunitcoeff` is an argument: `!iszero` over a field, `isunit` over
##  `Zn{n}`, and something narrower where only part of `Lambda` is inverted.
##
##  The search looks at the `limit` shortest basis elements.  Too small a
##  limit declines a genuine unit; it never accepts a non-unit, since the
##  candidate is verified.
##
"""
    monomialOf(h, isunitcoeff)

`(c, x)` where `h = c b_x` with `c` a unit, `nothing` otherwise.
"""
function monomialOf(h::PElt, isunitcoeff)
    (length(h.idx) == 1 && isunitcoeff(h.val[1])) || return nothing
    (h.val[1], Int(h.idx[1]))
end

##  b_x h = c b_y gives h^-1 = c^-1 b_y^-1 b_x, and h b_x = c b_y gives
##  h^-1 = c^-1 b_x b_y^-1
function unitSearch(A::PAlg{T}, h::PElt{T}, isunitcoeff, limit::Int) where T
    for k in 1:limit
        x = A.bylength[k]
        x == 1 && continue
        b = basis(A, x)
        m = monomialOf(mul(A, b, h), isunitcoeff)
        m === nothing ||
            return scale(inverse(m[1]), mul(A, basisInverse(A, m[2]), b))
        m = monomialOf(mul(A, h, b), isunitcoeff)
        m === nothing ||
            return scale(inverse(m[1]), mul(A, b, basisInverse(A, m[2])))
    end
    nothing
end

"""
    findUnitInverse(A, h, isunitcoeff; limit = 64)

`h^-1` where `h` is `c b`, `c b^-1 b'` or `c b b'^-1` for basis elements `b`,
`b'` and a unit `c` of `Lambda`; `nothing` otherwise, which does not say that
`h` is not a unit.  The inverse is verified on both sides before it is
returned.
"""
function findUnitInverse(A::PAlg{T}, h::PElt{T}, isunitcoeff;
                         limit::Integer = 64) where T
    iszero(h) && return nothing
    m = monomialOf(h, isunitcoeff)
    m === nothing || return scale(inverse(m[1]), basisInverse(A, m[2]))
    y = unitSearch(A, h, isunitcoeff, min(Int(limit), A.d))
    y === nothing ||
        (mul(A, h, y) == one(A) && mul(A, y, h) == one(A)) ||
        error("findUnitInverse: a wrong inverse for ", h)
    y
end

#############################################################################
##
##  The cyclic case.
##
##  These are the operations one generator has and several do not.  They
##  assert `iscyclic`, so that a general algebra handed to a caller that wants
##  `u^k` is an error and not a wrong answer.
##
function assertcyclic(A::PAlg, what)
    iscyclic(A) || throw(ArgumentError("$what: the algebra $(A.name) is not cyclic"))
end

"""
    ugen(A)

`u`, the generator of a cyclic `A` over `Lambda`.
"""
function ugen(A::PAlg{T}) where T
    assertcyclic(A, "ugen")
    A.d == 1 ? PElt{T}(Int32[1], T[A.params[1][1]]) : basis(A, 2)
end

"""
    uinv(A)

`u^-1`, from `c_0 u^-1 = u^(d-1) - c_(d-1) u^(d-2) - .. - c_1`.
"""
function uinv(A::PAlg{T}) where T
    assertcyclic(A, "uinv")
    v = mulinv(A, one(A), A.right[1], 1)
    mul(A, v, ugen(A)) == one(A) || error("u^-1 does not check out in $(A.name)")
    v
end

"""
    upow(A, k)

`u^k`, for `k` of either sign.
"""
function upow(A::PAlg{T}, k::Integer) where T
    assertcyclic(A, "upow")
    k >= 0 ? element(A, fill(1, k)) : rmulword(A, one(A), fill(-1, -k))
end

"""
    squareParams(params)

The minimal polynomial of `u^2` from that of `u`, `params` constant term
first, same degree `d`.

Splitting `m` into even and odd parts, `m(x) = E(x^2) + x O(x^2)`,

    prod_i (y - xi_i^2) = (-1)^d m(sqrt y) m(-sqrt y)
                        = (-1)^d (E(y)^2 - y O(y)^2)

which is the characteristic polynomial of multiplication by `u^2` and
annihilates `u^2` by Cayley-Hamilton.

The constant term is `(-1)^(d+1) c_0^2`, a unit exactly when `c_0` is.

It is the minimal polynomial unless two roots of `m` differ by a sign;
otherwise it is a proper multiple of it.
"""
function squareParams(params::AbstractVector{T}) where T
    d = length(params)
    d >= 1 || throw(ArgumentError("params is empty"))
    u = one(params[1])
    ##  a[k+1] is the coefficient of x^k in m(x) = x^d - sum c_k x^k
    a = [k == d ? u : -params[k+1] for k in 0:d]
    E = [a[k+1] for k in 0:d if iseven(k)]
    O = [a[k+1] for k in 0:d if isodd(k)]
    p = fill(zero(params[1]), d + 1)
    for i in eachindex(E), j in eachindex(E)
        p[i+j-1] += E[i] * E[j]
    end
    for i in eachindex(O), j in eachindex(O)
        p[i+j] -= O[i] * O[j]
    end
    isodd(d) && (p = -p)
    p[d+1] == u || error("the square polynomial is not monic: ", p[d+1])
    [-p[k+1] for k in 0:d-1]
end

end # module
