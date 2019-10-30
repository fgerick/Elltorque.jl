
mode2gradp(u,b,b0,Ω,ω) = coriolis(u,Ω) + lorentz(b,b0) - ω*u
mode2gradpmag(u,b,b0,Ω,ω) = coriolis(u,Ω) + Mire.advecterm(b,b0) + Mire.advecterm(b0,b) - ω*u
mode2gradpageo(u,ua,b,b0,Ω,ω) = coriolis(ua,Ω) + lorentz(b,b0) - ω*u

angularmom(u,cmatbulk,coordinate) = int_polynomial_ellipsoid((u × [x,y,z])[coordinate],cmatbulk)
emtorque(b,b0,cmatbulk,coordinate) = int_polynomial_ellipsoid(((Mire.advecterm(b,b0) + Mire.advecterm(b0,b)) × [x,y,z])[coordinate],cmatbulk)
emtorquejxb(b,b0,cmatbulk,coordinate) = int_polynomial_ellipsoid((Mire.lorentz(b,b0) × [x,y,z])[coordinate],cmatbulk)
hydropressuretorque(u,b,ω,b0,Ω,cmatbulk, coordinate) = int_polynomial_ellipsoid((mode2gradp(u,b,b0,Ω,ω) × [x,y,z])[coordinate],cmatbulk)
hydroageopressuretorque(u,b,ω,ua,b0,Ω,cmatbulk, coordinate) = int_polynomial_ellipsoid((mode2gradpageo(u,ua,b,b0,Ω,ω) × [x,y,z])[coordinate],cmatbulk)
totalpressuretorque(u,b,ω,b0,Ω,cmatbulk, coordinate) = int_polynomial_ellipsoid((mode2gradpmag(u,b,b0,Ω,ω) × [x,y,z])[coordinate],cmatbulk)
coriolistorque(u,Ω,cmatbulk,coordinate) = int_polynomial_ellipsoid((coriolis(u,Ω) × r)[coordinate],cmatbulk)


function torquebalance(N,a,b,c,vs,us,ug,ua,bs,ω,b0,Ω)
    cmatbulk = Mire.cacheint_Hsxyzpoly(N,a,b,c)
    cmatbulk3var = Mire.cacheint(N,a,b,c)*pi

    modes = [Mode(u,b,ω) for (u,b,ω) in zip(us,bs,ω)]
    cmatsurf = [cacheint_surface_torque(N+3,i,a,b,c) for i=1:3]

    Γp = [[hydropressuretorque(m,b0,Ω,cmatbulk3var,i) for m in modes] for i=1:3]
    Γpageo = [[hydroageopressuretorque(modes[j],ua[j],b0,Ω,cmatbulk,i) for j in 1:length(ua)] for i=1:3]

    Γptot = [[totalpressuretorque(m,b0,Ω,cmatbulk3var,i) for m in modes] for i=1:3]

    Γcor = [[coriolistorque(m.u,Ω,cmatbulk3var,i) for m in modes] for i=1:3]

    psmag = [sum(m.b.*b0) for m in modes]
    Γpmag = [[int_polynomial_ellipsoid(p,cmatsurf[i]) for p in psmag] for i=1:3]

    Lω = [[angularmom(m.u,cmatbulk3var,i) for m in modes].*ω for i=1:3]
    Lωa = [[angularmom(ua[j],cmatbulk,i) for j=1:length(ua)].*ω for i=1:3]

    Γem = [[emtorque(m.b,b0,cmatbulk3var,i) for m in modes] for i=1:3]

    return Γp, Γptot, Γpmag, Lω, Lωa, Γem, Γcor, Γpageo
end


# function loadandcalculatetorque(m::ModelSetup)
