using Revise, Mire, PyPlot, GenericSchur, DoubleFloats, JLD2,
      TypedPolynomials, MultivariatePolynomials, Printf, Arpack
using Destruct
CALCULATE = false
VERBOSE = true
SAVEFIGS = true
datapath = "../data/"
figpath = "../figs/"

PyPlot.rc("text",usetex=true)

# cd("/Users/gerickf/Repositories/pressure_torque_ellipsoid/run")
using Elltorque
cd(joinpath(dirname(pathof(Elltorque)),"../run"))
# include("../src/Elltorque.jl")


include("convergence.jl")
include("torques.jl")
# include("dense_angularmomentum.jl")
include("dense_sphere3d.jl")
include("ellipticity.jl")
include("lehnert.jl")
include("streamplots.jl")
include("malkus.jl")


#only seems to work in REPL mode...
# using Luxor
# include("ellipse_scetch.jl")
