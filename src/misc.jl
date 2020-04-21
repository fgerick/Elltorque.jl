# function to truncate the cartesian polynomials
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


#QG inertial mode frequencies following
# Zhang et al., J. Fluid Mech. (2001), vol. 437, pp. 103–119. eq (4.5)
ωZhangRossby(m::Int,N::Int)=-2/(m+2)*(√(1+m*(m+2)/(N*(2N+2m+1)))-1)+0.0im

# fast and slow mode frequencies following
# Malkus J. Fluid Mech. (1967), vol. 28, pp. 793-802, eq. (2.28)
ωMalkusSlow(m::Int, Le::Float64, N::Int, λ::T = ωZhangRossby(m,N)) where T<: Number = im*λ/2Le*(1 - √(1+4Le^2*m*(m-λ)/λ^2))
ωMalkusFast(m::Int, Le::Float64, N::Int, λ::T = ωZhangRossby(m,N)) where T<: Number = im*λ/2Le*(1 + √(1+4Le^2*m*(m-λ)/λ^2))
