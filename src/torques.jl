
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

    nm = length(us)
    Γp = [zeros(ComplexF64,nm) for i=1:3]
    Γpageo,Γptot,Γcor,Γpmag,Lω,Lωa,Γem = [deepcopy(Γp) for i=1:7]

    for i = 1:3
        cmatsurf = cacheint_surface_torque(N+3,i,a,b,c)
        for j = 1:nm
            Γp[i][j] = hydropressuretorque(us[j],bs[j],ω[j],b0,Ω,cmatbulk3var,i)
            Γpageo[i][j] = hydroageopressuretorque(us[j],bs[j],ω[j],ua[j],b0,Ω,cmatbulk,i)
            Γptot[i][j] = totalpressuretorque(us[j],bs[j],ω[j],b0,Ω,cmatbulk3var,i)
            Γcor[i][j] = coriolistorque(us[j],Ω,cmatbulk3var,i)
            Γpmag[i][j] = int_polynomial_ellipsoid(sum(bs[j].*b0),cmatsurf)
            Lω[i][j] = angularmom(us[j],cmatbulk3var,i)*ω[j]
            Lωa[i][j] = angularmom(ua[j],cmatbulk,i)*ω[j]
            Γem[i][j] = emtorque(bs[j],b0,cmatbulk3var,i)
        end
    end

    return Γp, Γptot, Γpmag, Lω, Lωa, Γem, Γcor, Γpageo
end

function loadandcalculatetorque(m::ModelSetup{T,D}, datapath="", SAVEDATA=false,dtypename="f64") where {T <: Real,D <: ModelDim}
    fname = joinpath(datapath,string(D)*"_$(m.name)_"*dtypename*"_N$(m.N).jld")
    JLD2.@load fname A B vs S ω evecs m Ω us bs
    ugua = split_ug_ua.(us,m.a,m.b,m.c)
    ug,ua = getindex.(ugua,1), getindex.(ugua,2)
    Γp, Γptot, Γpmag, Lω, Lωa, Γem, Γcor, Γpageo = torquebalance(m.N,m.a,m.b,m.c,vs,us,ug,ua,bs,ω,m.b0,Ω)
    if SAVEDATA
        fname = joinpath(datapath,"torquebalance_"*string(D)*"_$(m.name)_"*dtypename*"_N$(m.N).jld")
        JLD2.@save fname Γp Γptot Γpmag Lω Lωa Γem Γcor Γpageo
    end
    return true
end
