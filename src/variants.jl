#############################################################################
##
##  Variants of Relations  (Lists Version.)
##
module variants

export variantsRelations

function variantsRelations(genrel)
    abs(s::Int) = s < 0 ? genrel.invr[-s] : s
    variants = [Set{Vector{Int}}() for _ in genrel.gens]
    for (l, r) in genrel.rels
        relator =  vcat(l, -reverse(r))  # l * r^-1 = 1
        for (i, s) in enumerate(relator)
            vu = vcat(relator[i+1:end], relator[1:i-1])  #  u s v = 1
            push!(variants[abs(s)], abs.(-reverse(vu)))  #  s = (v u)^{-1}
            push!(variants[abs(-s)], abs.(vu))           #  s^-1 = v u
        end
    end
    return [sort(collect(v), by=length) for v in variants]
end

end # module
