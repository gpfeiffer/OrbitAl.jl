#############################################################################
##
#A  shifts.jl                                                         OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  Enumerate and visualize cyclic shifts
##
module shifts

export cyclic_shifts, cyclic_shifts_with_edges

using ..orbits
import ..coxeter.coxeterLength

function cyclic_shifts(W, w)
    function byCyclicShift(x, s)
        y = x^s
        coxeterLength(W, x) == coxeterLength(W, y) ?  y : x
    end
    orbit(W.gens, w, byCyclicShift)
end

function cyclic_shifts_with_edges(W, w)
    function byCyclicShift(x, s)
        y = x^s
        coxeterLength(W, x) == coxeterLength(W, y) ?  y : x
    end
    orbit_with_edges(W.gens, w, byCyclicShift)
end

end # module
