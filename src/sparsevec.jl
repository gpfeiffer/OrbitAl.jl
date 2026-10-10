#############################################################################
##
#A  sparsevec.jl                                                      OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  Sparse vectors: sorted positions and their nonzero coefficients
##
module sparsevec

export SparseVec, unitVec, drop

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
    SparseVec(i => c, ...)

The sparse vector with coefficient `c` at position `i`, for each pair.  The
coefficients of repeated positions are added up.
"""
function SparseVec(pairs::Pair...)
    poss = sort(unique(first.(pairs)))
    vals = [sum(c for (j, c) in pairs if j == i) for i in poss]
    return nonzero(poss, vals)
end

"""
    SparseVec(v::AbstractArray)

The sparse vector of the dense vector `v`.  A row vector, or any other array,
is read as the plain vector `v[:]`.
"""
function SparseVec(v::AbstractArray)
    v = v[:]
    poss = findall(!iszero, v)
    return SparseVec(poss, v[poss])
end

"""
    unitVec(K, i)

The unit vector ``e_i`` over `K`: a type, such as `Rational{Int}`, or a ring
that makes its own elements, such as a [`HeckeAlg`](@ref
OrbitAl.hecke.HeckeAlg).
"""
unitVec(K, i::Int) = SparseVec([i], [one(K)])

##  the type of the coefficients
Base.eltype(::Type{SparseVec{T}}) where T = T

Base.zero(::Type{SparseVec{T}}) where T = SparseVec(Int[], T[])
Base.zero(v::SparseVec) = zero(typeof(v))
Base.iszero(v::SparseVec) = isempty(v.poss)

##  the number of nonzero coefficients
Base.length(v::SparseVec) = length(v.poss)

##  the coefficient at position i, by binary search
function Base.getindex(v::SparseVec, i::Int)
    k = searchsortedfirst(v.poss, i)
    return k <= length(v.poss) && v.poss[k] == i ? v.vals[k] : zero(eltype(v))
end

"""
    drop(v::SparseVec, i)

The vector `v` without its coefficient at position `i`.
"""
function drop(v::SparseVec, i::Int)
    keep = v.poss .!= i
    return SparseVec(v.poss[keep], v.vals[keep])
end

##  merge the two sorted lists of positions, adding where they meet
function Base.:+(v::SparseVec, w::SparseVec)
    poss, vals = Int[], promote_type(eltype(v), eltype(w))[]
    j, k = 1, 1
    m, n = length(v), length(w)
    while j <= m || k <= n
        if k > n || j <= m && v.poss[j] < w.poss[k]        # only in v
            push!(poss, v.poss[j]); push!(vals, v.vals[j])
            j += 1
        elseif j > m || w.poss[k] < v.poss[j]               # only in w
            push!(poss, w.poss[k]); push!(vals, w.vals[k])
            k += 1
        else                                                # in both
            c = v.vals[j] + w.vals[k]
            iszero(c) || begin
                push!(poss, v.poss[j])
                push!(vals, c)
            end
            j += 1; k += 1
        end
    end
    return SparseVec(poss, vals)
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
    print(io, join([cf(c) * "e$i" for (i, c) in zip(v.poss, v.vals)], " + "))
end

end # module
