#functions to calculate the torques for a calulated mode

r = [x,y,z] #position vector

mode2gradp(u,b,b0,Ω,ω) = coriolis(u,Ω) + lorentz(b,b0) - ω*u
mode2gradpmag(u,b,b0,Ω,ω) = coriolis(u,Ω) + Mire.advecterm(b,b0) + Mire.advecterm(b0,b) - ω*u
mode2gradpageo(u,ua,b,b0,Ω,ω) = coriolis(ua,Ω) + lorentz(b,b0) - ω*u

angularmom(u,cmat,coordinate) = int_polynomial_ellipsoid((r × u)[coordinate],cmat)
emtorque(b,b0,cmat,coordinate) = int_polynomial_ellipsoid((r × (Mire.advecterm(b,b0) + Mire.advecterm(b0,b)))[coordinate],cmat)
hydropressuretorque(u,b,ω,b0,Ω,cmat, coordinate) = int_polynomial_ellipsoid((r × mode2gradp(u,b,b0,Ω,ω))[coordinate],cmat)
magpressuretorque(u,b,ω,b0,Ω,cmat, coordinate) = int_polynomial_ellipsoid((r × ∇(sum(b.*b0)))[coordinate],cmat)

totalpressuretorque(u,b,ω,b0,Ω,cmat, coordinate) = int_polynomial_ellipsoid((r × mode2gradpmag(u,b,b0,Ω,ω))[coordinate],cmat)
coriolistorque(u,Ω,cmat,coordinate) = int_polynomial_ellipsoid((r × coriolis(u,Ω))[coordinate],cmat)

#get all torques for all velocities us and associated magnetic field perturbations bs
#and other model setup properties (N,a,b,c,b0,Ω)
function torquebalance(N,a::T,b::T,c::T,us,bs,ω,b0,Ω; takeuni=false) where T
    cmat = Mire.cacheint(N,a,b,c)*pi
    nm = length(us)
    Γp = [zeros(Complex{T},nm) for i=1:3]
    Γpageo,Γptot,Γcor,Γpmag,Lω,Lωa,Γem = [deepcopy(Γp) for i=1:7]
    for i = 1:3
        for j = 1:nm
            if takeuni #projecting onto the uniform vorticity components
                u = projmode_uni(us[j],a,b,c,cmat, i)
            else
                u = us[j]
            end

            Γp[i][j] = hydropressuretorque(u,bs[j],ω[j],b0,Ω,cmat,i)
            Γptot[i][j] = totalpressuretorque(u,bs[j],ω[j],b0,Ω,cmat,i)
            Γcor[i][j] = coriolistorque(u,Ω,cmat,i)
            Γpmag[i][j] = magpressuretorque(u,bs[j],ω[j],b0,Ω,cmat,i)
            Lω[i][j] = angularmom(u,cmat,i)*ω[j]
            Γem[i][j] = emtorque(bs[j],b0,cmat,i)
        end
    end

    return Γp, Γptot, Γpmag, Lω, Γem, Γcor
end

# project a velocity u onto the i-th uniform vorticity component
function projmode_uni(u,a,b,c,cmat, i) where {T <: Real, D <: ModelDim}
    if i==1
        vx = [0,z/c^2,-y/b^2].* x^0*y^0*z^0
        vx /= √inner_product(cmat,vx,vx)
        v=vx
    elseif i==2
        vy = [z/c^2,0,-x/a^2].* x^0*y^0*z^0
        vy /= √inner_product(cmat,vy,vy)
        v=vy
    elseif i==3
        vz = [-y/b^2,x/a^2,0].* x^0*y^0*z^0
        vz /= √inner_product(cmat,vz,vz)
        v=vz
    else
        error("i=1,2,3")
    end
    return inner_product(cmat,u,v).*v
end


#convenience function to calculate the torque balance.
function loadandcalculatetorque(m::ModelSetup{T,D}, datapath="", SAVEDATA=false,dtypename="f64") where {T <: Real,D <: ModelDim}
    fname = joinpath(datapath,string(D)*"_$(m.name)_"*dtypename*"_N$(m.N).jld")
    JLD2.@load fname A B S ω evecs m Ω us bs

    Γp, Γptot, Γpmag, Lω, Γem, Γcor = torquebalance(m.N,m.a,m.b,m.c,us,bs,ω,m.b0,Ω)
    if SAVEDATA
        fname = joinpath(datapath,"torquebalance_"*string(D)*"_$(m.name)_"*dtypename*"_N$(m.N).jld")
        JLD2.@save fname Γp Γptot Γpmag Lω Γem Γcor
    end
    return true
end
