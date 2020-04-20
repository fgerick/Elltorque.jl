df641 = one(Double64)

a,b,c,Le = df64"1.25",df64"0.8",df641,df64"1e-5"

b0f = (a,b,c)->b0_1_1(a,b,c)+b0_1_3(a,b,c)
pAform = (x^0*y^0+x)/df64"3"
pAform2 = y^0*x/df64"3"

b0Af = (a,b,c)-> b0_Aform(pAform,a,b,c);

mqg = ModelSetup(a,b,c,Le,b0Af, "ellipse_aform7", 7, QG())
mhyb = ModelSetup(a,b,c,Le,b0_2_8, "ellipse_b289", 9, Hybrid());

if CALCULATE
    calculatemodes(mqg,datapath,true,"df64", ekin=false)
    calculatemodes(mhyb,datapath,true,"df64", ekin=false)
    loadandcalculatetorque(mqg,datapath,true,"df64")
    loadandcalculatetorque(mhyb,datapath,true,"df64")
end

function tbplot_z(m::ModelSetup{T,D},datapath="",figpath="",SAVEFIG=false;
                           dtypename = "df64",
                           ymax = 0,
                           xmax = 0,
                           lowlimexpx = -9,
                           lowlimexpy = -9,
                           ax = subplots(1,figsize=(4,2))[2],
                           symlogx = false,
                           lowlimlogx = m.Le
                           ) where {T <: Real, D <: ModelDim}
    fname = joinpath(datapath,"torquebalance_"*string(D)*"_$(m.name)_"*dtypename*"_N$(m.N).jld")
    JLD2.@load fname Γp Γptot Γpmag Lω Γem Γcor
    fname = joinpath(datapath,string(D)*"_$(m.name)_"*dtypename*"_N$(m.N).jld")
    JLD2.@load fname ω
    PyPlot.rc("text",usetex=true)

    finduniformmode = findmax(abs.(Lω[3]./ω))[2]

    onlypos = imag.(ω) .>= -eps()
    onlypos = (eachindex(ω).!=finduniformmode) .& onlypos

    ax.plot(abs.(imag.(ω[onlypos])), abs.(Γpmag[3][onlypos]) , linewidth=0, color="C2", marker="d", markersize=7,  label = L"|\Gamma_{\mathrm{pm},z}|") #"|∫b⋅b₀(n×r) dS|")
    ax.plot(abs.(imag.(ω[onlypos])), abs.(Γp[3][onlypos]) ,    linewidth=0, color="k", marker="+", markersize=8,  label = L"|\Gamma_{\mathrm{p},z}|") #"|∫(∇p×r) dV|")
    ax.plot(abs.(imag.(ω[onlypos])), abs.(Lω[3][onlypos]) ,    linewidth=0, color="C3", marker=".", markersize=10,  label = L"|\omega L_z|") #"|ω∫u×r dV|")
    ax.plot(abs.(imag.(ω[onlypos])), abs.(Γptot[3][onlypos]) , linewidth=0, color="k", marker="x", markersize=6,  label = L"|\Gamma_{\mathrm{p},z}+\Gamma_{\mathrm{pm},z}|")#"|∫∇(p+b⋅b₀)×r dV|")

    ax.plot([abs.(imag.(ω[finduniformmode]))], [abs.(Γpmag[3][finduniformmode])] , linewidth=0, color="C2", marker="d",alpha=0.25, markersize=7 ) #"|∫b⋅b₀(n×r) dS|")
    ax.plot([abs.(imag.(ω[finduniformmode]))], [abs.(Γp[3][finduniformmode])] ,    linewidth=0, color="k", marker="+",alpha=0.25, markersize=8 ) #"|∫(∇p×r) dV|")
    ax.plot([abs.(imag.(ω[finduniformmode]))], [abs.(Lω[3][finduniformmode])] ,    linewidth=0, color="C3", marker=".",alpha=0.25 ,markersize=10 ) #"|ω∫u×r dV|")
    ax.plot([abs.(imag.(ω[finduniformmode]))], [abs.(Γptot[3][finduniformmode])] , linewidth=0, color="k", marker="x",alpha=0.25, markersize=6  )#"|∫∇(p+b⋅b₀)×r dV|")

    ax.set_ylabel(L"|\Gamma_z|")
    ax.set_yscale("symlog", linthreshy=10.0^lowlimexpy)

    yl = ax.get_ylim()[2]
    ylm = ymax!=0 ? ymax : maximum(yl) #max(ymax,maximum(yl))

    ax.set_yticks(vcat([0],10.0.^(lowlimexpy:3:ceil(log10(ylm)))))
    ax.set_ylim([0,10.0^ceil(log10(ylm))])

    xl = ax.get_xlim()[2]

    xl = xmax!=0 ? xmax : xl #max(ymax,maximum(yl))
    if symlogx
        ax.set_xscale("symlog", linthreshx=10.0^lowlimexpx)
        ax.set_xticks(vcat([0],10.0.^(lowlimexpx:3:(ceil(log10(xl))))))
        ax.set_xlim([0,10.0^ceil(log10(xl))])
    else
        ax.set_xscale("log")
        ax.set_xlim([lowlimlogx,10.0^ceil(log10(xl))])

    end

    if SAVEFIG
        figname="torquebalance_z_"*string(D)*"_$(m.name)_"*dtypename*"_N$(m.N).pdf"
        f.savefig(joinpath(figpath,figname), bbox_inches="tight")
    end
end


#get numerical values of torsional modes torques:
function tbplot_z_tmvalues(m::ModelSetup{T,D},datapath="",figpath="",SAVEFIG=false;
                           dtypename = "df64") where {T <: Real, D <: ModelDim}
    fname = joinpath(datapath,"torquebalance_"*string(D)*"_$(m.name)_"*dtypename*"_N$(m.N).jld")
    JLD2.@load fname Γp Γptot Γpmag Lω Γem Γcor
    fname = joinpath(datapath,string(D)*"_$(m.name)_"*dtypename*"_N$(m.N).jld")
    JLD2.@load fname ω
    PyPlot.rc("text",usetex=true)

    tm = eachindex(ω)[0.1.<imag.(ω).<10]
    for i in tm
        @printf "Im(ω) = %.3f , ωL = %.3f + %.3fi, Γp = %.3f + %.3fi, Γpm = %.3f + %.3fi \n" imag(ω[i]) real(-Lω[3][i]) imag(-Lω[3][i]) real(Γp[3][i]) imag(Γp[3][i]) real(Γpmag[3][i]) imag(Γpmag[3][i])
        @printf "Γp + Γpm = %.2f + %.2fi, Γptot = %.2f + %.2fi \n" real(Γp[3][i]+Γpmag[3][i]) imag(Γp[3][i]+Γpmag[3][i]) real(Γptot[3][i]) imag(Γptot[3][i])
    end
end



#Figure 9 of the paper

f,ax=subplots(1,figsize=(3,2.2))
aval=0.4
cm=PyPlot.get_cmap(:cividis);
ax.fill_between(range(1e-1,10,length=10),0,100,color=cm(0.5),alpha=aval,linewidth=0,label=L"\mathrm{TM}")
ax.fill_between(range(1e-6,1e-2,length=10),0,100,color=cm(0),alpha=aval,linewidth=0,label=L"\mathrm{slow}")
ax.fill_between(range(1e2,1e6,length=10),0,100,color=cm(1.0),alpha=aval,linewidth=0,label=L"\mathrm{fast}")

tbplot_z(mqg, datapath, figpath, false, lowlimexpx = -6, lowlimexpy = -10,ymax=10,xmax=10^6,ax=ax,lowlimlogx=1e-6)
ax.legend(ncol=2,loc="center",bbox_to_anchor=(0.5,1.3))
ax.set_xlabel(L"|\omega|")

if SAVEFIGS
    figname="torquebalance_z_qg.pdf"
    f.savefig(joinpath(figpath,figname), bbox_inches="tight")
end




# Figure A2 of the paper

f,ax=subplots(1,figsize=(3,2.2))
aval=0.4
cm=PyPlot.get_cmap(:cividis);
ax.fill_between(range(1e-1,10,length=10),0,100,color=cm(0.5),alpha=aval,linewidth=0,label=L"\mathrm{TM}")
ax.fill_between(range(1e-6,1e-2,length=10),0,100,color=cm(0),alpha=aval,linewidth=0,label=L"\mathrm{slow}")
ax.fill_between(range(1e2,1e6,length=10),0,100,color=cm(1.0),alpha=aval,linewidth=0,label=L"\mathrm{fast}")
tbplot_z(mhyb, datapath, figpath, false, lowlimexpx = -6, lowlimexpy = -10,ymax=10,xmax=10^6,ax=ax,lowlimlogx=1e-6)
ax.legend(ncol=2,loc="center",bbox_to_anchor=(0.5,1.3))
ax.set_xlabel(L"|\omega|")
if SAVEFIGS
    figname="torquebalance_z_hybrid.pdf"
    f.savefig(joinpath(figpath,figname), bbox_inches="tight")
end


#Numerical values for Earth like ellipticities:


ϵ = df64"1e-3"

b = ((1 - ϵ)/(1 + ϵ))^(1//4)
a = 1/b

mqg_ϵ = ModelSetup(a,b,c,df64"1e-5",b0Af, "ellipse_aform7_eps", 7, QG())
mhyb_ϵ = ModelSetup(a,b,c,df64"1e-5",b0_2_8, "ellipse_b289_eps", 9, Hybrid());

if CALCULATE
    calculatemodes(mqg_ϵ,datapath,true,"df64", ekin=true)
    calculatemodes(mhyb_ϵ,datapath,true,"df64", ekin=true)
    loadandcalculatetorque(mqg_ϵ,datapath,true,"df64")
    loadandcalculatetorque(mhyb_ϵ,datapath,true,"df64")
end


ϵ = df64"1e-4"
b = ((1 - ϵ)/(1 + ϵ))^(1//4)
a = 1/b
mqg_ϵ4 = ModelSetup(a,b,c,df64"1e-5",b0Af, "ellipse_aform7_eps4", 7, QG())
mhyb_ϵ4 = ModelSetup(a,b,c,df64"1e-5",b0_2_8, "ellipse_b289_eps4", 9, Hybrid());


if CALCULATE
    calculatemodes(mqg_ϵ4,datapath,true,"df64", ekin=true)
    calculatemodes(mhyb_ϵ4,datapath,true,"df64", ekin=true)
    loadandcalculatetorque(mqg_ϵ4,datapath,true,"df64")
    loadandcalculatetorque(mhyb_ϵ4,datapath,true,"df64")
end


using Latexify, Unitful


function tm_values_1m(m0::ModelSetup{T,D},u,iTM,datapath="";
                           dtypename = "df64") where {T <: Real, D <: ModelDim}

    fname = joinpath(datapath,"torquebalance_"*string(D)*"_$(m0.name)_"*dtypename*"_N$(m0.N).jld")
    JLD2.@load fname Γp Γptot Γpmag Lω Γem Γcor
    fname = joinpath(datapath,string(D)*"_$(m0.name)_"*dtypename*"_N$(m0.N).jld")
    JLD2.@load fname ω


    tm_inds = eachindex(ω)[1e-1.<imag.(ω).<10]

    A = Unitful.A
    kg = Unitful.kg
    s = Unitful.s
    m = Unitful.m
    K = Unitful.K
    yr = 365.25*24*3600s
    Ω = kg*m^2/(s^3*A^2)
    H = Ω*s
    Tesla = kg/(A*s^2)
    W = kg*m^2/(s^3)
    Pa = kg/(m*s^2)
    km = 1e3*m
    μ0 = Unitful.μ0

    Ωrot = 2pi/(24*3600*s)
    ρ    = 1.2e4*kg/m^3
    r0   = 3478e3*m

    ω1  = abs.(ω[tm_inds][iTM])
    ω6y = 2π/6yr


    Le = ω6y/(ω1*Ωrot)
    B0 =upreferred(Le*Ωrot*r0*√(ρ*μ0))

    u0 = u*m/s
    vA = upreferred(B0/sqrt(μ0*ρ))
    U = u0/sqrt(3/8pi)
    Nm = kg*m^2/s^2

    torque_fact = r0^4*U*ρ*2pi/6yr
    L = Lω[3][tm_inds[iTM]]/ω[tm_inds[iTM]]
    ωLz = abs(Lω[3][tm_inds[iTM]]/ω[tm_inds[iTM]]*ω6y*U*ρ*r0^4/Nm)
    ω   = imag(ω[tm_inds[iTM]])
    return ellipticity(m0),upreferred(1e3*B0/(Tesla)),Le,ω,ωLz,abs(L)
end



dim(m::ModelSetup{T,D}) where {T <: Real, D <: ModelDim} = D
ellipticity(m::ModelSetup{T,D}) where {T <: Real, D <: ModelDim} = (m.a^2-m.b^2)/(m.a^2+m.b^2)

function tm_values_all(ms,u,iTM,datapath="";
                           dtypename = "df64")

    tab = Array{Any}(undef,length(ms)*length(iTM)+1,6)
    tab[1,:] = [L"\mathrm{Model}",L"\omega [2\pi/T_A]",L"L_z",L"\mathrm{Le}",L"B_0\, [\mathrm{mT}]",L"\omega L_z\, [\mathrm{Nm}]"]
    for (i,m) in enumerate(ms)
        for (j,itm) in enumerate(iTM)
            ϵ,B0,Le,ω,ωLz,L = tm_values_1m(m,u,itm,datapath,dtypename = dtypename)
            tab[((i-1)*length(iTM)+1)+j,:] = [LaTeXString("\\mathrm{"*string(dim(m))*"}"),ω,L*10^5,Le,B0,ωLz]
        end
    end
   latexify(tab,env=:table, fmt=FancyNumberFormatter(3, "\\times"))
end


tm_values_all([mqg_ϵ,mhyb_ϵ],5e-6,3:-1:1,datapath)
