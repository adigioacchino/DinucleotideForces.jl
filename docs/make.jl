using DinucleotideForces
using Documenter

DocMeta.setdocmeta!(DinucleotideForces, :DocTestSetup, :(using DinucleotideForces); recursive=true)

makedocs(;
    modules=[DinucleotideForces],
    authors="Andrea Di Gioacchino",
    repo="https://github.com/adigioacchino/DinucleotideForces.jl/blob/{commit}{path}#{line}",
    sitename="DinucleotideForces.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://adigioacchino.github.io/DinucleotideForces.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/adigioacchino/DinucleotideForces.jl",
    devbranch="main",
)
