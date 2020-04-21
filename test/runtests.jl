using Test, Elltorque

global CALCULATE = true
global VERBOSE = true
global SAVEFIGS = true
global ellpath = joinpath(dirname(pathof(Elltorque)),"..")
global runpath = joinpath(ellpath,"run")
global datapath = joinpath(ellpath,"data/")
global figpath = joinpath(ellpath,"figs/")
try
    mkdir(datapath)
    mkdir(joinpath(datapath,"runlehnert_b0An11/"))
catch
end

const GROUP = get(ENV, "GROUP","all")

PyPlot.rc("text",usetex=true)
PyPlot.pygui(false)

if GROUP == "convergence" || GROUP == "all"
    @testset "convergence" begin; @test include("../run/convergence.jl"); end

elseif GROUP == "lehnert" || GROUP == "all"
    @testset "lehnert" begin; @test include("../run/lehnert.jl"); end

elseif GROUP == "ellipticity" || GROUP == "all"
    @testset "ellipticity" begin; @test include("../run/ellipticity.jl"); end

elseif GROUP == "torques" || GROUP == "all"
    @testset "torques" begin
        @test include("../run/torques.jl")
        @test include("../run/streamplots.jl")
    end

elseif GROUP == "dense_amom" || GROUP == "all"
    @testset "dense angular momentum" begin; @test include("../run/dense_angularmomentum.jl"); end

elseif GROUP == "dense_sphere" || GROUP == "all"
    @testset "dense sphere" begin; @test include("../run/dense_sphere3d.jl"); end
    
elseif GROUP == "malkus" || GROUP == "all"
    @testset "malkus" begin; @test include("../run/malkus.jl"); end
end

# @testset "Calculate and plot" begin
#     @test Elltorque.run(true)
# end
