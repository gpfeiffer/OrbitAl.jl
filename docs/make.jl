using Documenter
using OrbitAl
using OrbitAl.coxeter
using OrbitAl.simsgroup
using OrbitAl.syt
using OrbitAl.coset
using OrbitAl.involution

makedocs(
    sitename = "OrbitAl.jl",
    modules = [OrbitAl, OrbitAl.coxeter, OrbitAl.simsgroup, OrbitAl.syt, OrbitAl.coset, OrbitAl.involution],
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
        ],
    ],
    authors = "Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>",
)

deploydocs(
    repo = "github.com/gpfeiffer/OrbitAl.jl.git",
    devbranch = "main",
)
