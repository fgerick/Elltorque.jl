using Pkg
Pkg.add(PackageSpec(url="https://github.com/fgerick/Mire.jl.git"))
Pkg.develop(PackageSpec(path=joinpath(@__DIR__)))
