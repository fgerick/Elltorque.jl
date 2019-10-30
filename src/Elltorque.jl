module Elltorque

using Reexport, LinearAlgebra, Statistics, GenericSchur, DoubleFloats, JLD2
@reexport using Mire

export ModelSetup, ModelDim, Full, Hybrid

abstract type ModelDim; end

struct Full <: ModelDim; end
struct Hybrid <: ModelDim; end

# struct QG <: ModelDim; end TODO


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

export b0_1_1,b0_1_2,b0_1_3,b0_2_6,b0_3_6

include("analysis.jl")
export calculatemodes

include("geosplit.jl")
export split_ug_ua

include("torques.jl")
export torquebalance

end # module
