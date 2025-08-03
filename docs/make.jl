using Documenter
using OrbitAl

# DocMeta.setdocmeta!(OrbitAl, :DocTestSetup, :(using OrbitAl); recursive=true)

makedocs(
    sitename = "OrbitAl.jl",
    modules = [OrbitAl],
    checkdocs = :none,
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "Permutations" => "permutation.md",
        "Permutation Groups" => "permgroup.md",
        "Orbits" => "orbits.md",
        "Coxeter Groups" => "coxeter.md",
    ],
    authors = "Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>",
)
