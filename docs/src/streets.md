# Streets and Descent Algebras

!!! note
    This is an on-demand module. Load it with `using OrbitAl.streets`.

The alleys of a finite Coxeter group ``W`` with simple reflections ``S`` are the
paths in the Hasse diagram of the power set of ``S``; the streets are the
orbits of a conjugation action of the free monoid ``S^*`` on the alleys.  The
streets span the street algebra ``\Xi``, a path algebra, which maps onto
Solomon's descent algebra ``\Sigma(W)`` by the anti-homomorphism ``\Delta``.
See G. Pfeiffer, *A quiver presentation for Solomon's descent algebra*,
Adv. Math. 220 (2009) 1428--1465.

Everything is combinatorial: from ``W`` only the opposition involutions
``s \mapsto w_J s w_J`` of the parabolic subsets ``J`` are used.  Subsets are
bit masks, and the alleys are numbered without a table, so that the rank can
be as large as memory allows (type ``D_9``, with ``2\,681\,216`` alleys, takes
about ten seconds).

```@docs
OrbitAl.streets.StreetSystem
OrbitAl.streets.shapeReps
OrbitAl.streets.xiQuiver
OrbitAl.streets.descentQuiver
```
