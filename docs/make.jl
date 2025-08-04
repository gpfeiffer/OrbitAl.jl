using Documenter
using OrbitAl

makedocs(
    sitename = "OrbitAl.jl",
    modules = [OrbitAl],
    checkdocs = :none,
    doctest = true,
    format = Documenter.HTML(),
    repo = Remotes.GitHub("gpfeiffer", "OrbitAl.jl"),
    pages = [
        "Home" => "index.md",
        "Permutations" => "permutation.md",
        "Orbits" => "orbits.md",
    ],
    authors = "Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>",
)

deploydocs(
    repo = "github.com/gpfeiffer/OrbitAl.jl.git",
    devbranch = "main",
)
