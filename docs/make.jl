using Documenter, Elltorque

makedocs(;
    modules=[Elltorque],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/fgerick/Elltorque.jl/blob/{commit}{path}#L{line}",
    sitename="Elltorque.jl",
    authors="Felix Gerick",
    assets=String[],
)

deploydocs(;
    repo="github.com/fgerick/Elltorque.jl",
)
