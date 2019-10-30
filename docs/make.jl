using Documenter, Elltorque

makedocs(;
    modules=[Elltorque],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    sitename="Elltorque.jl",
    authors="Felix Gerick",
    assets=String[],
)

deploydocs(;
    repo="github.com/fgerick/Elltorque.jl",
)
