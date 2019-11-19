function truncpoly(u;atol=√(eps()))
    uo=zero(u)
    ti = terms(u)
    for t in ti
        c=coefficient(t)
        if abs(c)>atol
            uo+=t
        end
    end
    return uo
end

function truncvec(u;atol=√(eps()))
    uo=zero(u)
    for i=1:3
        ti = terms(u[i])
        for t in ti
            c=coefficient(t)
            if abs(c)>atol
                uo[i]+=t
            end
        end
    end
    return uo
end



ωZhangRossby(m::Int,N::Int)=-2/(m+2)*(√(1+m*(m+2)/(N*(2N+2m+1)))-1)+0.0im


ωMalkusSlow(m::Int, Le::Float64, N::Int, λ::T = ωZhangRossby(m,N)) where T<: Number = im*λ/2Le*(1 - √(1+4Le^2*m*(m-λ)/λ^2))

"""
Fast inertial wave frequencies λ following Malkus (1967). (See Labbe et al. 2015, eq. 23)
"""
ωMalkusFast(m::Int, Le::Float64, N::Int, λ::T = ωZhangRossby(m,N)) where T<: Number = im*λ/2Le*(1 + √(1+4Le^2*m*(m-λ)/λ^2))
