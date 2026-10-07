#############################################################################
##
#A  sparsevec.jl                                                      OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  Sparse vectors: sorted positions and their nonzero coefficients
##
module sparsevec

export SparseVec, unitVec

"""
    SparseVec{T}

A sparse vector over `T`: the positions `poss` of its nonzero coefficients,
in increasing order, and the coefficients `vals`.  A `SparseVec` has no
length and its positions can be any positive integers.
"""
struct SparseVec{T}
    poss::Vector{Int}  # the positions of the nonzero coefficients, increasing
    vals::Vector{T}    # the coefficients
end

##  drop the zero coefficients
nonzero(poss, vals) = (nz = .!iszero.(vals); SparseVec(poss[nz], vals[nz]))

"""
    SparseVec(p => c, ...)

The sparse vector with coefficient `c` at position `p`, for each pair.  The
coefficients of repeated positions are added up.
"""
function SparseVec(pairs::Pair...)
    poss = sort(unique(first.(pairs)))
    vals = [sum(c for (q, c) in pairs if q == p) for p in poss]
    return nonzero(poss, vals)
end

"""
    unitVec(T, p)

The unit vector ``e_p`` over `T`.
"""
unitVec(T::Type, p::Int) = SparseVec([p], [one(T)])

Base.zero(::Type{SparseVec{T}}) where T = SparseVec(Int[], T[])
Base.zero(v::SparseVec) = zero(typeof(v))
Base.iszero(v::SparseVec) = isempty(v.poss)

##  the number of nonzero coefficients
Base.length(v::SparseVec) = length(v.poss)

##  the coefficient at position p, by binary search
function Base.getindex(v::SparseVec, p::Int)
    k = searchsortedfirst(v.poss, p)
    return k <= length(v.poss) && v.poss[k] == p ? v.vals[k] : zero(eltype(v.vals))
end

function Base.:+(v::SparseVec, w::SparseVec)
    poss = sort(union(v.poss, w.poss))
    return nonzero(poss, [v[p] + w[p] for p in poss])
end

Base.:*(c, v::SparseVec) = nonzero(v.poss, c .* v.vals)
Base.:-(v::SparseVec) = SparseVec(v.poss, -v.vals)
Base.:-(v::SparseVec, w::SparseVec) = v + (-w)

##  division by a unit c
Base.:/(v::SparseVec, c) = SparseVec(v.poss, v.vals ./ c)

Base.:(==)(v::SparseVec, w::SparseVec) = v.poss == w.poss && v.vals == w.vals
Base.hash(v::SparseVec, h::UInt) = hash(v.poss, hash(v.vals, h))

function Base.show(io::IO, v::SparseVec)
    length(v) > 0 || return print(io, "0")
    cf(c) = isone(c) ? "" : "(" * sprint(show, c; context = :typeinfo => typeof(c)) * ")"
    print(io, join([cf(c) * "e$p" for (p, c) in zip(v.poss, v.vals)], " + "))
end

end # module
