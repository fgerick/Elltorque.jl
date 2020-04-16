
function integral_u1(i,j,k,a,b,c)
    if iseven(i) && isodd(i+j) && iseven(k)
        return -4*s^(i+j-1)*H^k*a^(i-1)*b^j*c^k*gamma((1 + i)/2)*gamma(1 + j/2)/(2*(1+k)*gamma((3 + i + j)/2))
    else
        return zero(a)
    end
end

function integral_u2(i,j,k,a,b,c)
    if isodd(i) && isodd(i+j) && iseven(k)
        return 4*s^(i+j-1)*H^k*a^i*b^(j-1)*c^k*gamma(1 + i/2)*gamma((1 + j)/2)/(2*(1+k)*gamma((3 + i + j)/2))
    else
        return zero(a)
    end
end

function int_polynomial(p,intfun::Function,a,b,c)
    cs=coefficients(p)
    m=monomials(p)
    e=exponents.(m)
    out=0
    for i=1:length(cs)
        out+=cs[i]*intfun(e[i]...,a,b,c)
    end
    return out
end

u_g(u,a,b,c) = (int_polynomial(u[1],integral_u1,a,b,c) + int_polynomial(u[2],integral_u2,a,b,c))/2π .* [-y*a/b,x*b/a,0]

function split_ug_ua(u,a,b,c)
    mvarhack = H^0*s^0*x^0*y^0*z^0
    ug = u_g(u,a,b,c) .*mvarhack
    ua = (u - ug).*mvarhack
    return ug, ua
end
