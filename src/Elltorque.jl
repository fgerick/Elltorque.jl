module Elltorque

using Reexport
@reexport using Mire, LinearAlgebra, Statistics, GenericSchur, DoubleFloats, JLD2,
      MultivariatePolynomials, TypedPolynomials, Distributed, LinearMaps,
      Arpack, PyPlot, Destruct, Printf, Latexify, Unitful, SparseArrays




abstract type ModelDim; end

struct Full <: ModelDim; end
struct Hybrid <: ModelDim; end
struct QG <: ModelDim; end

struct ModelSetup{T <: Real, D <: ModelDim}
    a::T
    b::T
    c::T
    Le::T
    b0
    name::String
    N::Int
end



ModelSetup(a::T,b::T,c::T,Le::T,b0fun::Function,name::String,N::Int, ::D) where {T <: Real, D <: ModelDim} = ModelSetup{T,D}(a,b,c,Le,b0fun(a,b,c),name,N)

#selected mean magnetic fields satisfying b₀⋅n=0 and ∇⋅b₀=0 from Wu & Roberts, 2011.


#linear base
b0_1_1(a,b,c) = [0,-z/c^2,y/b^2]
b0_1_2(a,b,c) = [z/c^2,0,-x/a^2]
b0_1_3(a,b,c) = [-y/b^2,x/a^2,0]

#quadratic base
b0_2_6(a,b,c) = [-z*x,z*y,c^2*(x^2/a^2-y^2/b^2)]

#cubic base
b0_3_6(a,b,c) = [-z*x^2,2*z*x*y,c^2*x*(x^2/a^2-2y^2/b^2)]

# b₀⋅n = 0 at δV only:

fijkl(i,j,k,l,a,b,c) = i - j*x^2/a^2 - k*y^2/b^2 - l*z^2/c^2

b0_2_7(a,b,c) = [a^2*fijkl(1,1,2,2,a,b,c), x*y, x*z]
b0_2_8(a,b,c) = [y*x, b^2*fijkl(1,2,1,2,a,b,c), z*y]
b0_2_9(a,b,c) = [z*x, z*y, c^2*fijkl(1,2,2,1,a,b,c)]

#Aformulation

function b0_Aform(p_xy,a,b,c)
    ts = terms(p_xy)
    out = zero(Mire.qg_vel(0,0,a,b,c))
    for t in ts
        out .+= coefficient(t).*Mire.qg_vel(exponents(t)...,a,b,c)
    end
    return out
end


include("misc.jl")

include("analysis.jl")

include("torques.jl")

include("track_malkus.jl")

"""
    run(calculate=false)

Produces all plots except the ellipsoid scetch (Figure 2). If calculate is set
to `true` all data will be computed.

"""
function run(calculate=false)
    if calculate
        global CALCULATE = true
    else
        global CALCULATE = false
    end
    global VERBOSE = true
    global SAVEFIGS = true
    global ellpath = joinpath(dirname(pathof(Elltorque)),"..")
    global runpath = joinpath(ellpath,"run")
    global datapath = joinpath(ellpath,"data/")
    try
        mkdir(datapath)
        mkdir(joinpath(datapath,"runlehnert_b0An11/"))
    catch

    end
    global figpath = joinpath(ellpath,"figs/")

    PyPlot.rc("text",usetex=true)

    PyPlot.pygui(false)
    include(joinpath(runpath,"convergence.jl"))
    println("convergence done.")
    include(joinpath(runpath,"torques.jl"))
    println("torques done.")
    include(joinpath(runpath,"dense_angularmomentum.jl"))
    println("dense angular momentum done.")
    include(joinpath(runpath,"dense_sphere3d.jl"))
    println("dense sphere done.")
    include(joinpath(runpath,"ellipticity.jl"))
    println("ellipticity done.")
    close("all")
    include(joinpath(runpath,"lehnert.jl"))
    println("lehnert done.")
    include(joinpath(runpath,"streamplots.jl"))
    println("streamplots done.")
    include(joinpath(runpath,"malkus.jl"))
    println("malkus done.")
    close("all")
    println("done.")
    return true
end

function plot_scetch()
    cd(joinpath(dirname(pathof(Elltorque)),"../run"))
    include("ellipse_scetch.jl")
end


export ModelSetup, ModelDim, Full, Hybrid, QG, calculatemodes, split_ug_ua,
       torquebalance, loadandcalculatetorque, runcalculations,tracking_lehnert,
       tracking_ellipt, tracking_ellipt_reverse, get_ub, get_ub1, track_malkus

export angularmom, mode2gradp, emtorque, hydropressuretorque, magpressuretorque,
        totalpressuretorque, coriolistorque
export b0_1_1,b0_1_2,b0_1_3,b0_2_6,b0_2_7,b0_2_8,b0_2_9,b0_3_6, b0_Aform
#
end # module
