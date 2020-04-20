using Pkg
Pkg.add(PackageSpec(url="https://github.com/fgerick/Mire.jl.git"))
Pkg.dev(PackageSpec(path=joinpath(@__DIR__)))
