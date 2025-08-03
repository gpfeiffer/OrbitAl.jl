using Documenter
using OrbitAl

makedocs(
    sitename = "OrbitAl.jl",
    modules = [OrbitAl],
    format = Documenter.HTML(),
    repo = "github.com/gpfeiffer/OrbitAl.jl",
    devbranch = "main",
    pages = [
        "Home" => "index.md",
        "Permutations" => "permutations.md",
        "Orbits" => "orbits.md",
    ],
    authors = "Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>",
)

deploydocs(
    repo = "github.com/gpfeiffer/OrbitAl.jl.git",
    devbranch = "main",
)
