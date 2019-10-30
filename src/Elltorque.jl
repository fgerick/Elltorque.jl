module Elltorque

using LinearAlgebra, Statistics, Mire, GenericSchur, DoubleFloats

struct ModelSetup{T <: Real}
    a::T
    b::T
    c::T
    Le::T
    b0
end

ModelSetup(a::T,b::T,c::T,Le::T,b0fun::Function) where T <: Real = ModelSetup{T}(a,b,c,Le,b0fun(a,b,c))

#selected mean magnetic fields satisfying b₀⋅n=0 and ∇⋅b₀=0 from Wu & Roberts, 2011.




end # module
