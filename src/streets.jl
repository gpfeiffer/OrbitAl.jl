#############################################################################
##
#A  streets.jl                                                        OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  Alleys, streets, the street algebra Xi and the quiver of the descent
#C  algebra, after G. Pfeiffer, A quiver presentation for Solomon's descent
#C  algebra, Adv. Math. 220 (2009) 1428--1465 [Pf09].
##
##  Everything here is combinatorial: the only input taken from the Coxeter
##  group W is the opposition involution  s -> w_J s w_J  of each parabolic
##  subset J.  The conjugation of J by omega(J, s) = w_J w_L, L = J + s, is
##  then  x -> opp_L(opp_J(x)).
##
##  Subsets of S = 1:n are bit masks (bit s-1 for s).  An alley (L; s_1, ...,
##  s_l) is a source L and a sequence of distinct elements of L; the alleys
##  are numbered 1, 2, ... by L, then by l, then by the mixed radix number of
##  the positions of the s_i among the remaining elements of L.  So there is
##  no table of alleys: `alley_id` and `decode!` translate.
##
module streets

using ..permutation
using ..coxeter: longestElt

export StreetSystem, nrAlleys, nrStreets, nrShapes, shapeReps
export xiQuiver, descentQuiver, descentPresentation

#############################################################################
##
##  The street system: alley numbering, streets, shapes.
##
struct StreetSystem
    n::Int                       # the rank, |S|
    opp::Vector{Vector{Int8}}    # opp[J+1][s] = w_J s w_J, for s in J
    off::Vector{Vector{Int}}     # off[L+1][l+1]: alleys before (L, length l)
    nalleys::Int
    street::Vector{Int32}        # the street of each alley
    rep::Vector{Int}             # the smallest alley of each street
    src::Vector{Int32}           # the source shape of each street
    tgt::Vector{Int32}           # the target shape of each street
    len::Vector{Int8}            # the length of each street
    shapeOf::Vector{Int32}       # shapeOf[J+1]: the shape of the subset J
    shapeRep::Vector{Int}        # a subset (mask) in each shape
    starts::Vector{Int}          # the offsets of the (L, l) ranges, sorted,
    startL::Vector{Int}          # with their sources
    startl::Vector{Int}          # and lengths
end

bit(s) = 1 << (s - 1)
elements(J, n) = [s for s in 1:n if J & bit(s) != 0]

##  the opposition involution of each parabolic subset, from W acting on roots
function oppositions(W)
    n = W.data[:rank];  N = W.data[:N]
    opp = Vector{Vector{Int8}}(undef, 1 << n)
    for J in 0:(1 << n) - 1
        o = zeros(Int8, n)
        if J != 0
            wJ = longestElt(W, elements(J, n))
            for s in 1:n
                J & bit(s) != 0 && (o[s] = s^wJ - N)    # w_J a_s = -a_{opp(s)}
            end
        end
        opp[J + 1] = o
    end
    return opp
end

##  k!/(k-l)!: the number of alleys of length l with a source of size k
arrangements(k, l) = prod(k-l+1:k; init = 1)

##  the number of the alley (L; seq[1:l])
function alley_id(off, L, seq, l)
    k = count_ones(L);  rem = L;  r = 0
    for i in 1:l
        s = seq[i]
        r = r * (k - i + 1) + count_ones(rem & (bit(s) - 1))
        rem &= ~bit(s)
    end
    return off[L + 1][l + 1] + r + 1
end

##  the i-th (from 0) element of the subset rem
function nth_element(rem, d)
    for s in 1:64
        if rem & bit(s) != 0
            d == 0 && return s
            d -= 1
        end
    end
    error("no such element")
end

##  write the segments of the alley with number r (from 0) among those of
##  source L and length l into buf
function decode!(buf, L, l, r)
    k = count_ones(L)
    for i in l:-1:1
        q, d = divrem(r, k - i + 1)
        buf[i] = d;  r = q
    end
    rem = L
    for i in 1:l
        s = nth_element(rem, buf[i])
        buf[i] = s;  rem &= ~bit(s)
    end
    return buf
end

##  the source, length and rank (from 0) of the alley with number id
function locate(X, id)
    i = searchsortedlast(X.starts, id - 1)
    return X.startL[i], X.startl[i], id - 1 - X.starts[i]
end

##  iterate over all alleys: f(id, L, l, seq)
function foreach_alley(f, n, off)
    buf = zeros(Int, n)
    for L in 0:(1 << n) - 1
        k = count_ones(L)
        for l in 0:k
            o = off[L + 1][l + 1]
            for r in 0:arrangements(k, l) - 1
                decode!(buf, L, l, r)
                f(o + r + 1, L, l, buf)
            end
        end
    end
end

##  the image of a subset under the conjugation x -> opp_M(opp_L(x)), L ⊆ M
function conj_set(opp, L, M, J, n)
    oL = opp[L + 1];  oM = opp[M + 1];  K = 0
    for s in 1:n
        J & bit(s) != 0 && (K |= bit(oM[oL[s]]))
    end
    return K
end

function uf_find(parent, i)
    r = i
    while parent[r] != r
        r = parent[r]
    end
    while parent[i] != r        # path compression
        parent[i], i = r, parent[i]
    end
    return r
end

"""
    StreetSystem(W)

The alleys, streets and shapes of the finite Coxeter group `W`
(a `CoxeterGp`).  The streets are the orbits of the conjugation action of the
free monoid on the alleys [Pf09, Section 6]; the shapes are the streets of
length 0.

# Examples
```jldoctest
julia> using OrbitAl, OrbitAl.coxeter, OrbitAl.streets

julia> X = StreetSystem(CoxeterGp(cartanMat("A", 3)));

julia> (nrAlleys(X), nrStreets(X), nrShapes(X))
(38, 27, 5)
```
"""
function StreetSystem(W)
    n = W.data[:rank]
    n <= 16 || throw(ArgumentError("rank $n is too large"))
    opp = oppositions(W)
    off = Vector{Vector{Int}}(undef, 1 << n);  total = 0
    for L in 0:(1 << n) - 1
        k = count_ones(L);  o = Int[]
        for l in 0:k
            push!(o, total);  total += arrangements(k, l)
        end
        off[L + 1] = o
    end
    total < typemax(Int32) || throw(ArgumentError("too many alleys: $total"))

    # the streets: Union-Find along the edges a -> a.r of the action graph
    parent = collect(Int32, 1:total)
    img = zeros(Int, n)
    foreach_alley(n, off) do id, L, l, seq
        oL = opp[L + 1]
        for r in 1:n
            L & bit(r) != 0 && continue
            oM = opp[(L | bit(r)) + 1]
            for i in 1:l
                img[i] = oM[oL[seq[i]]]
            end
            jd = alley_id(off, conj_set(opp, L, L | bit(r), L, n), img, l)
            a, b = uf_find(parent, Int32(id)), uf_find(parent, Int32(jd))
            a == b && continue
            a < b ? (parent[b] = a) : (parent[a] = b)
        end
    end
    street = zeros(Int32, total);  rep = Int[]
    for id in 1:total
        r = uf_find(parent, Int32(id))
        if r == id
            push!(rep, id);  street[id] = length(rep)
        else
            street[id] = street[r]
        end
    end
    parent = nothing

    # the shapes, numbered in order of their smallest subset
    shapeOf = zeros(Int32, 1 << n);  shapeRep = Int[];  shapeNr = Dict{Int32,Int32}()
    for J in 0:(1 << n) - 1
        st = street[off[J + 1][1] + 1]
        shapeOf[J + 1] = get!(shapeNr, st) do
            push!(shapeRep, J);  length(shapeRep)
        end
    end
    m = length(rep)
    src = zeros(Int32, m);  tgt = zeros(Int32, m);  len = zeros(Int8, m)
    foreach_alley(n, off) do id, L, l, seq
        g = street[id]
        if rep[g] == id
            T = L
            for i in 1:l
                T &= ~bit(seq[i])
            end
            src[g] = shapeOf[L + 1];  tgt[g] = shapeOf[T + 1];  len[g] = l
        end
    end
    starts = Int[];  startL = Int[];  startl = Int[]
    for L in 0:(1 << n) - 1, l in 0:count_ones(L)
        push!(starts, off[L + 1][l + 1]);  push!(startL, L);  push!(startl, l)
    end
    return StreetSystem(n, opp, off, total, street, rep, src, tgt, len,
                        shapeOf, shapeRep, starts, startL, startl)
end

nrAlleys(X::StreetSystem) = X.nalleys
nrStreets(X::StreetSystem) = length(X.rep)
nrShapes(X::StreetSystem) = length(X.shapeRep)

"""
    shapeReps(X)

A representative subset of each shape of `X`, as a sorted list of
simple reflections.
"""
shapeReps(X::StreetSystem) = [elements(J, X.n) for J in X.shapeRep]

Base.show(io::IO, X::StreetSystem) = print(io, "StreetSystem of rank $(X.n): ",
    "$(nrAlleys(X)) alleys, $(nrStreets(X)) streets, $(nrShapes(X)) shapes")

#############################################################################
##
##  Linear algebra mod p.
##
const P = 2147483647    # 2^31 - 1

##  the rank of a list of sparse rows (Dict position => coefficient) mod P
function rank_mod_p(rows)
    piv = Dict{Int, Dict{Int,Int}}()
    for v0 in rows
        v = Dict{Int,Int}(k => mod(x, P) for (k, x) in v0 if mod(x, P) != 0)
        while !isempty(v)
            c = minimum(keys(v))
            if haskey(piv, c)
                f = v[c]
                for (k, x) in piv[c]
                    y = mod(get(v, k, 0) - f * x, P)
                    y == 0 ? delete!(v, k) : (v[k] = y)
                end
            else
                inv = powermod(v[c], P - 2, P)
                piv[c] = Dict(k => mod(x * inv, P) for (k, x) in v)
                break
            end
        end
    end
    return length(piv)
end

#############################################################################
##
##  The street algebra Xi and its radical.
##
##  The product of streets is a 0/1 combination of streets: gamma occurs in
##  alpha o beta iff the prefix of length l(alpha) of one (hence each) alley
##  in gamma lies in alpha and the remaining suffix in beta.  So rad^2 Xi is
##  spanned by the sets  {gamma : split at k gives (alpha, beta)}.
##
##  the streets of each block (lam, mu), lam != mu
function blocks(X::StreetSystem)
    B = Dict{Tuple{Int32,Int32}, Vector{Int32}}()
    for g in 1:nrStreets(X)
        X.len[g] > 0 && push!(get!(() -> Int32[], B, (X.src[g], X.tgt[g])), g)
    end
    return B
end

##  the generators of rad^2 Xi in one block, as sets of streets
function rad2_generators(X::StreetSystem, gs)
    n = X.n
    gens = Dict{NTuple{3,Int32}, Vector{Int32}}()
    seq = zeros(Int, n);  suf = zeros(Int, n)
    for g in gs
        L, l, r = locate(X, X.rep[g])
        decode!(seq, L, l, r)
        K = L
        for k in 1:l-1
            K &= ~bit(seq[k])
            p = alley_id(X.off, L, seq, k)
            for i in k+1:l
                suf[i - k] = seq[i]
            end
            q = alley_id(X.off, K, suf, l - k)
            push!(get!(() -> Int32[], gens, (Int32(k), X.street[p], X.street[q])), g)
        end
    end
    return values(gens)
end

##  the number of paths from lam to mu in the quiver with arrow counts A
function path_counts(X::StreetSystem, A)
    m = nrShapes(X)
    order = sortperm([-count_ones(J) for J in X.shapeRep])
    cnt = zeros(Int, m, m)
    for lam in 1:m
        cnt[lam, lam] = 1
        for mu in order
            mu == lam && continue
            cnt[lam, mu] = sum(cnt[lam, nu] * get(A, (nu, mu), 0) for nu in 1:m)
        end
    end
    return cnt
end

"""
    xiQuiver(X)

The Ext-quiver of the street algebra ``\\Xi`` of `X`: a named tuple with the
number of arrows `arrows[(lam, mu)]` from shape `lam` to shape `mu` (in the
direction of the streets, from source to target), the Cartan invariants
`cartan[(lam, mu)]`, and `ispath`, which says whether ``\\Xi`` is the path
algebra of this quiver, i.e., whether the Cartan invariants count its paths.
(It is in all cases computed so far; that it always is, i.e., that ``\\Xi``
is hereditary, is [Pf09, Conjecture 6.14].)

# Examples
```jldoctest
julia> using OrbitAl, OrbitAl.coxeter, OrbitAl.streets

julia> q = xiQuiver(StreetSystem(CoxeterGp(cartanMat("A", 3))));

julia> (sum(values(q.arrows)), q.ispath)
(8, true)
```
"""
function xiQuiver(X::StreetSystem)
    cartan = Dict{Tuple{Int32,Int32},Int}()
    for g in 1:nrStreets(X)
        key = (X.src[g], X.tgt[g])
        cartan[key] = get(cartan, key, 0) + 1
    end
    arrows = Dict{Tuple{Int32,Int32},Int}()
    for (key, gs) in blocks(X)
        c = length(gs)
        r2 = rank_mod_p([Dict(Int(x) => 1 for x in s) for s in rad2_generators(X, gs)])
        c > r2 && (arrows[key] = c - r2)
    end
    cnt = path_counts(X, arrows)
    m = nrShapes(X)
    ispath = all(cnt[lam, mu] == get(cartan, (lam, mu), 0) for lam in 1:m, mu in 1:m)
    return (arrows = arrows, cartan = cartan, ispath = ispath)
end

#############################################################################
##
##  The descent algebra, via the anti-homomorphism Delta: Xi -> Sigma(W).
##
##  Delta(a) = Delta(b) - Delta(b)^omega, for b = sigma(a) and omega =
##  omega(K_s, s), where K is the source of a and s its first segment
##  [Pf09, Section 7 and Proposition 7.2].  So Delta(a) is computed along the
##  suffixes of a, without a table.  Only alleys whose source is the chosen
##  representative L of its shape are used: f_alpha e_L = f_{L o alpha}
##  [Pf09, Corollary 9.3], and x -> x e_L is injective on Sigma e_lambda.
##
function delta(X::StreetSystem, L, l, seq)
    n = X.n
    K = L
    for i in 1:l
        K &= ~bit(seq[i])
    end
    v = Dict{Int,Int}(K => 1)        # Delta of the target, a path of length 0
    Ks = K
    for k in l:-1:1
        Kk = Ks | bit(seq[k])        # the source of the suffix starting at k
        w = Dict{Int,Int}()
        for (J, x) in v
            w[J] = get(w, J, 0) + x
            J2 = conj_set(X.opp, Ks, Kk, J, n)
            w[J2] = get(w, J2, 0) - x
        end
        v = filter(p -> p.second != 0, w)
        Ks = Kk
    end
    return v
end

##  Delta(L o gamma) for each street gamma, with L the representative of its
##  source, as a sparse vector over the subsets; the alleys with source L are
##  consecutive
function delta_reps(X::StreetSystem)
    n = X.n
    D = [Pair{Int32,Int}[] for _ in 1:nrStreets(X)]
    seq = zeros(Int, n)
    for L in X.shapeRep
        acc = Dict{Int32, Dict{Int,Int}}()
        for l in 1:count_ones(L)
            o = X.off[L + 1][l + 1]
            for r in 0:arrangements(count_ones(L), l) - 1
                decode!(seq, L, l, r)
                a = get!(() -> Dict{Int,Int}(), acc, X.street[o + r + 1])
                for (J, x) in delta(X, L, l, seq)
                    a[J] = get(a, J, 0) + x
                end
            end
        end
        for (g, a) in acc
            D[g] = [Int32(J) => x for (J, x) in a if x != 0]
        end
    end
    return D
end

todict(v) = Dict{Int,Int}(Int(J) => x for (J, x) in v)

function addup(D, gs)
    v = Dict{Int,Int}()
    for g in gs, (J, x) in D[g]
        v[J] = get(v, J, 0) + x
    end
    return v
end

"""
    descentQuiver(X)

The quiver of the descent algebra ``\\Sigma(W)`` of `X`, computed with
[Pf09, Algorithm 10.1]: a named tuple with the number of arrows
`arrows[(lam, mu)]` for each pair of shapes (in the direction of the
streets, from source `lam` to target `mu`; [Pf09] draws them from `mu` to
`lam`), and `dimrad`, the dimension of the radical, which is
``2^n - |\\Lambda|``.

# Examples
```jldoctest
julia> using OrbitAl, OrbitAl.coxeter, OrbitAl.streets

julia> q = descentQuiver(StreetSystem(CoxeterGp(cartanMat("A", 3))));

julia> (sum(values(q.arrows)), q.dimrad)
(2, 3)
```
"""
function descentQuiver(X::StreetSystem)
    D = delta_reps(X)
    arrows = Dict{Tuple{Int32,Int32},Int}();  dimrad = 0
    for (key, gs) in blocks(X)
        r1 = rank_mod_p([todict(D[g]) for g in gs])
        r2 = rank_mod_p([addup(D, s) for s in rad2_generators(X, gs)])
        dimrad += r1
        r1 > r2 && (arrows[key] = r1 - r2)
    end
    return (arrows = arrows, dimrad = dimrad)
end

#############################################################################
##
##  Incremental elimination mod P, with the combinations of the inputs.
##
mutable struct Eliminator
    piv::Dict{Int, Tuple{Dict{Int,Int}, Dict{Int,Int}}}   # pivot => (row, comb)
end
Eliminator() = Eliminator(Dict{Int, Tuple{Dict{Int,Int}, Dict{Int,Int}}}())

axpy!(v, f, w) = for (k, x) in w
    y = mod(get(v, k, 0) + f * x, P)
    y == 0 ? delete!(v, k) : (v[k] = y)
end

##  reduce v (with combination c) by the pivots; add it as a new pivot unless
##  it reduces to 0.  Returns the reduced combination if v reduces to 0, and
##  nothing otherwise.
function eliminate!(E::Eliminator, v0, c0)
    v = Dict{Int,Int}(k => mod(x, P) for (k, x) in v0 if mod(x, P) != 0)
    c = Dict{Int,Int}(k => mod(x, P) for (k, x) in c0 if mod(x, P) != 0)
    while !isempty(v)
        k = minimum(keys(v))
        if haskey(E.piv, k)
            f = P - v[k];  (w, d) = E.piv[k]
            axpy!(v, f, w);  axpy!(c, f, d)
        else
            inv = powermod(v[k], P - 2, P)
            E.piv[k] = (Dict(i => mod(x * inv, P) for (i, x) in v),
                        Dict(i => mod(x * inv, P) for (i, x) in c))
            return nothing
        end
    end
    return c
end

##  a/b with |a|, b <= sqrt(P/2) and a/b = x mod P
function ratrecon(x)
    x = mod(x, P);  x == 0 && return 0//1
    r0, r1, s0, s1 = P, x, 0, 1
    bound = isqrt(P ÷ 2)
    while r1 > bound
        q = r0 ÷ r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    end
    abs(s1) <= bound || return nothing
    return (sign(s1) * r1) // abs(s1)
end

##  a basis of the span of the vectors v (Dicts over 1:N) in reduced echelon
##  form, pivots at the largest positions; coefficients as rationals
function echelon_rational(vs)
    rows = Dict{Int, Dict{Int,Int}}()
    for v0 in vs
        v = Dict{Int,Int}(k => mod(x, P) for (k, x) in v0 if mod(x, P) != 0)
        while !isempty(v)
            k = maximum(keys(v))
            haskey(rows, k) || break
            axpy!(v, P - v[k], rows[k])
        end
        isempty(v) && continue
        k = maximum(keys(v));  inv = powermod(v[k], P - 2, P)
        for i in keys(v)
            v[i] = mod(v[i] * inv, P)
        end
        for (j, w) in rows                      # keep it reduced
            haskey(w, k) && axpy!(w, P - w[k], v)
        end
        rows[k] = v
    end
    return [Dict(i => ratrecon(x) for (i, x) in rows[k]) for k in sort(collect(keys(rows)))]
end

##  Delta(J o x) for all J in the shape of the source L of x, from v =
##  Delta(L o x): conjugate along the action graph of the shape, using
##  Delta(a^y) = Delta(a)^y [Pf09, Proposition 7.2] and J o x = (L o x)^y
function spread(X::StreetSystem, L, v)
    rows = Dict{Int, Dict{Int,Int}}(L => v);  todo = [L]
    while !isempty(todo)
        K = pop!(todo)
        for s in 1:X.n
            K & bit(s) != 0 && continue
            M = K | bit(s);  K2 = conj_set(X.opp, K, M, K, X.n)
            haskey(rows, K2) && continue
            rows[K2] = Dict(conj_set(X.opp, K, M, J, X.n) => x for (J, x) in rows[K])
            push!(todo, K2)
        end
    end
    return rows
end

##  the reduced echelon form of the vectors (with their combinations), pivots
##  at the smallest positions, sorted by pivot
function rref_with_combs(E::Eliminator)
    ks = sort(collect(keys(E.piv)))
    rows = Dict(k => (copy(E.piv[k][1]), copy(E.piv[k][2])) for k in ks)
    for k in reverse(ks), k2 in ks
        k2 == k && continue
        (w, d) = rows[k2]
        f = get(w, k, 0)
        f == 0 && continue
        axpy!(w, P - f, rows[k][1]);  axpy!(d, P - f, rows[k][2])
    end
    return [rows[k] for k in ks]
end

"""
    descentPresentation(X; arrows = :orthogonal)

A presentation of the descent algebra ``\\Sigma(W)`` of `X` as a quiver with
relations, by [Pf09, Algorithm 10.1].

The arrows from `lam` to `mu` are elements of ``\\Xi`` whose images span a
complement of ``\\Delta(\\mathrm{rad}^2 \\Xi)`` in ``\\Delta(\\mathrm{rad}\\, \\Xi)``
in that block.  With `arrows = :orthogonal` (the default), this is the
orthogonal complement for the inner product in which the basis ``e_J`` of
``\\Sigma(W)`` [Pf09, Proposition 9.2] is orthonormal (restricted to
``\\Sigma(W) e_L`` for the representative ``L`` of `lam`), with the basis in
reduced echelon form in the coordinates ``e_J``.  This choice does not
depend on a choice of streets.  With `arrows = :streets`, each arrow is the
first street of its block that is independent modulo
``\\Delta(\\mathrm{rad}^2 \\Xi)``.

Returns a named tuple with

* `arrows`: a list of triples `(comb, source, target)`, one for each arrow of
  the quiver, in the direction of the streets (from source to target), where
  `comb` is the arrow as a linear combination of streets (a `Dict`, with
  coefficients mod ``2^{31} - 1``);
* `paths`: for each pair `(lam, mu)` of shapes, the list of all paths from
  `lam` to `mu`, as lists of arrow numbers;
* `relations`: for each pair `(lam, mu)`, a basis of the relations among the
  paths from `lam` to `mu` that are not consequences of relations between
  shorter paths, in reduced echelon form, as `Dict`s from path numbers to
  rational coefficients (`nothing` if rational reconstruction fails);
* `allrelations`: the dimensions of the full spaces of relations.

Paths compose like streets: a path `[e1, e2]` is the product
``\\alpha_1 \\circ \\alpha_2`` in ``\\Xi``, which ``\\Delta`` maps to
``f_{\\alpha_2} f_{\\alpha_1}``, as ``\\Delta`` is an anti-homomorphism.

# Examples
```jldoctest
julia> using OrbitAl, OrbitAl.coxeter, OrbitAl.streets

julia> P = descentPresentation(StreetSystem(CoxeterGp(cartanMat("A", 5))));

julia> (length(P.arrows), sum(length, values(P.relations)))
(9, 1)
```
"""
function descentPresentation(X::StreetSystem; arrows = :orthogonal)
    arrows in (:orthogonal, :streets) ||
        throw(ArgumentError("arrows must be :orthogonal or :streets"))
    mode = arrows
    D = delta_reps(X)

    # the arrows, with their images Delta(L o arrow)
    arrows = Tuple{Dict{Int32,Int},Int32,Int32}[];  avec = Dict{Int,Int}[]
    for (key, gs) in sort(collect(blocks(X)))
        U = Eliminator()
        for v in (addup(D, s) for s in rad2_generators(X, gs))
            eliminate!(U, v, Dict{Int,Int}())
        end
        if mode == :streets
            for g in gs
                v = todict(D[g])
                if isnothing(eliminate!(U, v, Dict{Int,Int}()))
                    push!(arrows, (Dict{Int32,Int}(g => 1), key[1], key[2]))
                    push!(avec, Dict(k => mod(x, P) for (k, x) in v if mod(x, P) != 0))
                end
            end
            continue
        end
        V = Eliminator()                      # a basis of Delta(block), with streets
        for g in gs
            eliminate!(V, todict(D[g]), Dict(Int(g) => 1))
        end
        bs = collect(values(V.piv))
        us = [w for (w, _) in values(U.piv)]
        length(bs) == length(us) && continue
        # the combinations of the bs orthogonal to all us
        G = Eliminator()
        perp = Dict{Int,Int}[]
        for (i, (b, _)) in enumerate(bs)
            g = Dict{Int,Int}()
            for (j, u) in enumerate(us)
                x = 0
                for (k, y) in b
                    x = mod(x + y * get(u, k, 0), P)
                end
                x == 0 || (g[j] = x)
            end
            c = eliminate!(G, g, Dict(i => 1))
            isnothing(c) || push!(perp, c)
        end
        length(perp) == length(bs) - length(us) ||
            error("degenerate inner product mod P in block $key")
        C = Eliminator()
        for c in perp
            w = Dict{Int,Int}();  comb = Dict{Int,Int}()
            for (i, x) in c
                axpy!(w, x, bs[i][1]);  axpy!(comb, x, bs[i][2])
            end
            eliminate!(C, w, comb)
        end
        for (w, comb) in rref_with_combs(C)
            push!(arrows, (Dict{Int32,Int}(Int32(g) => x for (g, x) in comb), key[1], key[2]))
            push!(avec, w)
        end
    end

    # mu(arrow): the subset J maps to Delta(J o arrow)
    mu = [spread(X, X.shapeRep[a[2]], avec[i]) for (i, a) in enumerate(arrows)]
    function act(v, i)
        w = Dict{Int,Int}()
        for (J, x) in v, (K, y) in get(mu[i], J, Dict{Int,Int}())
            w[K] = mod(get(w, K, 0) + x * y, P)
        end
        return filter(p -> p.second != 0, w)
    end

    # all paths, by depth-first search, with their images
    out = Dict{Int32, Vector{Int}}()
    for (i, a) in enumerate(arrows)
        push!(get!(() -> Int[], out, a[2]), i)
    end
    paths = Dict{Tuple{Int32,Int32}, Vector{Vector{Int}}}()
    images = Dict{Tuple{Int32,Int32}, Vector{Dict{Int,Int}}}()
    function dfs(lam, path, v)
        mu_ = arrows[path[end]][3]
        push!(get!(() -> Vector{Int}[], paths, (lam, mu_)), copy(path))
        push!(get!(() -> Dict{Int,Int}[], images, (lam, mu_)), v)
        for j in get(out, mu_, Int[])
            push!(path, j);  dfs(lam, path, act(v, j));  pop!(path)
        end
    end
    for (i, a) in enumerate(arrows)
        dfs(a[2], [i], avec[i])
    end

    # the relations in each block, and the consequences of shorter ones
    kernel = Dict{Tuple{Int32,Int32}, Vector{Dict{Int,Int}}}()
    index = Dict(key => Dict(p => i for (i, p) in enumerate(ps)) for (key, ps) in paths)
    for (key, ps) in paths
        E = Eliminator();  ker = Dict{Int,Int}[]
        for (i, v) in enumerate(images[key])
            c = eliminate!(E, v, Dict(i => 1))
            isnothing(c) || push!(ker, c)
        end
        kernel[key] = ker
    end
    relations = Dict{Tuple{Int32,Int32}, Vector{Dict{Int,Rational{Int}}}}()
    allrelations = Dict{Tuple{Int32,Int32}, Int}()
    for (key, ker) in kernel
        isempty(ker) && continue
        lam, mu_ = key;  idx = index[key]
        cons = Dict{Int,Int}[]
        for (i, a) in enumerate(arrows)          # e . rho and rho . e
            if a[2] == lam
                for rho in get(kernel, (a[3], mu_), Dict{Int,Int}[])
                    ps = paths[(a[3], mu_)]
                    push!(cons, Dict(idx[[i; ps[j]]] => x for (j, x) in rho))
                end
            end
            if a[3] == mu_
                for rho in get(kernel, (lam, a[2]), Dict{Int,Int}[])
                    ps = paths[(lam, a[2])]
                    push!(cons, Dict(idx[[ps[j]; i]] => x for (j, x) in rho))
                end
            end
        end
        E = Eliminator()
        for c in cons
            eliminate!(E, c, Dict{Int,Int}())
        end
        new = Dict{Int,Int}[]
        for rho in ker
            isnothing(eliminate!(E, rho, Dict{Int,Int}())) && push!(new, rho)
        end
        allrelations[key] = length(ker)
        isempty(new) || (relations[key] = echelon_rational(new))
    end
    return (arrows = arrows, paths = paths, relations = relations,
            allrelations = allrelations)
end

end # module streets
