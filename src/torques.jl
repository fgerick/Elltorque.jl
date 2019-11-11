
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


function torquebalance(N,a::T,b::T,c::T,us,ug,ua,bs,ω,b0,Ω; takeuni=false) where T
    cmatbulk = real.(Mire.cacheint_Hsxyzpoly(N,a,b,c))
    cmatbulk3var = Mire.cacheint(N,a,b,c)*pi
    nm = length(us)
    Γp = [zeros(Complex{T},nm) for i=1:3]
    Γpageo,Γptot,Γcor,Γpmag,Lω,Lωa,Γem = [deepcopy(Γp) for i=1:7]
    for i = 1:3
        cmatsurf = cacheint_surface_torque(N+3,i,a,b,c)
        for j = 1:nm
            if takeuni
                u = projmode_uni(us[j],a,b,c,cmatbulk3var, i)
            else
                u = us[j]
            end

            Γp[i][j] = hydropressuretorque(u,bs[j],ω[j],b0,Ω,cmatbulk3var,i)
            # Γpageo[i][j] = hydroageopressuretorque(us[j],bs[j],ω[j],ua[j],b0,Ω,cmatbulk,i)
            Γptot[i][j] = totalpressuretorque(u,bs[j],ω[j],b0,Ω,cmatbulk3var,i)
            Γcor[i][j] = coriolistorque(u,Ω,cmatbulk3var,i)
            Γpmag[i][j] = int_polynomial_ellipsoid(sum(bs[j].*b0),cmatsurf)
            Lω[i][j] = angularmom(u,cmatbulk3var,i)*ω[j]
            # Lωa[i][j] = angularmom(ua[j],cmatbulk,i)*ω[j]
            Γem[i][j] = emtorque(bs[j],b0,cmatbulk3var,i)
        end
    end

    # return Γp, Γptot, Γpmag, Lω, Lωa, Γem, Γcor, Γpageo
    return Γp, Γptot, Γpmag, Lω, Γem, Γcor
end

function projmode_uni(u,a,b,c,cmatbulk, i) where {T <: Real, D <: ModelDim}
    if i==1
        vx = [0,z/c^2,-y/b^2].* x^0*y^0*z^0
        vx /= √inner_product(cmatbulk,vx,vx)
        v=vx
    elseif i==2
        vy = [z/c^2,0,-x/a^2].* x^0*y^0*z^0
        vy /= √inner_product(cmatbulk,vy,vy)
        v=vy
    elseif i==3
        vz = [-y/b^2,x/a^2,0].* x^0*y^0*z^0
        vz /= √inner_product(cmatbulk,vz,vz)
        v=vz
    else
        error("i=1,2,3")
    end
    return inner_product(cmatbulk,u,v).*v
end

function univortangmom(m::ModelSetup{T,D},datapath=""; dtypename = "df64") where {T <: Real, D <: ModelDim}

    fname = joinpath(datapath,string(D)*"_$(m.name)_"*dtypename*"_N$(m.N).jld")
        if D == Full
            JLD2.@load fname A B vs S ω evecs m Ω us bs
        elseif D == Hybrid
            JLD2.@load fname A B vs vs_qg S ω evecs m Ω us bs
        elseif D == QG
            JLD2.@load fname A B vs_qg S ω evecs m Ω us bs
        end

    nev = length(ω)
    a,b,c = m.a,m.b,m.c
    cmatbulk = cacheint(m.N,a,b,c)*pi

    vx = [0,z/c^2,-y/b^2].* x^0*y^0*z^0
    vy = [z/c^2,0,-x/a^2].* x^0*y^0*z^0
    vz = [-y/b^2,x/a^2,0].* x^0*y^0*z^0

    vx /= √inner_product(cmatbulk,vx,vx)
    vy /= √inner_product(cmatbulk,vy,vy)
    vz /= √inner_product(cmatbulk,vz,vz)


    if D == Full
        vs_uni = [vx,vy,vz]
    else
        vs_uni = [zero(vz),zero(vz),vz]
    end

    us_uni = [[] for i=1:3]
    L = [zeros(Complex{T},nev) for i=1:3]
    for i = 1:3
        v_uni = vs_uni[i]
        for j = 1:nev
            u_uni = inner_product(cmatbulk,us[j],v_uni).*v_uni
            L[i][j] = int_polynomial_ellipsoid((u_uni × [x,y,z])[i],cmatbulk)
            push!(us_uni[i],u_uni)
        end
    end
    return L, us_uni
end

#
# function loadandcalculatetorque(m::ModelSetup{T,D}, datapath="", SAVEDATA=false,dtypename="f64") where {T <: Real,D <: ModelDim}
#     # @everywhere begin
#     fname = joinpath(datapath,string(D)*"_$(m.name)_"*dtypename*"_N$(m.N).jld")
#     JLD2.@load fname A B S ω evecs m Ω us bs
#
#     ugua = split_ug_ua.(us,m.a,m.b,m.c)
#     ug,ua = getindex.(ugua,1), getindex.(ugua,2)
#     # end
#     Γp, Γptot, Γpmag, Lω, Lωa, Γem, Γcor, Γpageo = torquebalance(m.N,m.a,m.b,m.c,us,ug,ua,bs,ω,m.b0,Ω)
#     if SAVEDATA
#         fname = joinpath(datapath,"torquebalance_"*string(D)*"_$(m.name)_"*dtypename*"_N$(m.N).jld")
#         JLD2.@save fname Γp Γptot Γpmag Lω Lωa Γem Γcor Γpageo
#     end
#     return true
# end
#
# function loadandcalculatetorquethread(mutex, m::ModelSetup{T,D}, datapath="", SAVEDATA=false,dtypename="f64") where {T <: Real,D <: ModelDim}
#     # @everywhere begin
#     fname = joinpath(datapath,string(D)*"_$(m.name)_"*dtypename*"_N$(m.N).jld")
#     JLD2.@load fname A B S ω evecs m Ω us bs
#
#     ugua = split_ug_ua.(us,m.a,m.b,m.c)
#     ug,ua = getindex.(ugua,1), getindex.(ugua,2)
#     # end
#     Γp, Γptot, Γpmag, Lω, Lωa, Γem, Γcor, Γpageo = torquebalance(m.N,m.a,m.b,m.c,us,ug,ua,bs,ω,m.b0,Ω)
#     lock(mutex)
#     if SAVEDATA
#         fname = joinpath(datapath,"torquebalance_"*string(D)*"_$(m.name)_"*dtypename*"_N$(m.N).jld")
#         JLD2.@save fname Γp Γptot Γpmag Lω Lωa Γem Γcor Γpageo
#     end
#     unlock(mutex)
#     return true
# end
