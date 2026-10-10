using Documenter
using OrbitAl
using OrbitAl.coxeter
using OrbitAl.simsgroup
using OrbitAl.syt
using OrbitAl.coset
using OrbitAl.involution
using OrbitAl.modn
using OrbitAl.sparsevec
using OrbitAl.unionfind
using OrbitAl.hecke
using OrbitAl.linear
using OrbitAl.enumerator
using OrbitAl.vectorenum

makedocs(
    sitename = "OrbitAl.jl",
    modules = [OrbitAl, OrbitAl.coxeter, OrbitAl.simsgroup, OrbitAl.syt, OrbitAl.coset, OrbitAl.involution,
               OrbitAl.modn, OrbitAl.sparsevec, OrbitAl.unionfind, OrbitAl.hecke, OrbitAl.linear, OrbitAl.enumerator, OrbitAl.vectorenum],
    checkdocs = :none,
    doctest = true,
    format = Documenter.HTML(),
    repo = Remotes.GitHub("gpfeiffer", "OrbitAl.jl"),
    pages = [
        "Home" => "index.md",
        "Permutations" => "permutation.md",
        "Orbits" => "orbits.md",
        "Permutation Groups" => "permgroup.md",
        "On Demand" => [
            "Coxeter Groups" => "coxeter.md",
            "Schreier-Sims Groups" => "simsgroup.md",
            "Standard Young Tableaux" => "syt.md",
            "Cosets" => "coset.md",
            "Involutions" => "involution.md",
            "Modular Arithmetic" => "modn.md",
            "Sparse Vectors" => "sparsevec.md",
            "Union-Find" => "unionfind.md",
            "Hecke Algebras" => "hecke.md",
            "Linear Orbit Algorithms" => "linear.md",
            "Coset Enumeration" => "enumerator.md",
            "Vector Enumeration" => "vectorenum.md",
            "Example Groups" => "data.md",
        ],
    ],
    authors = "Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>",
)

deploydocs(
    repo = "github.com/gpfeiffer/OrbitAl.jl.git",
    devbranch = "main",
)
