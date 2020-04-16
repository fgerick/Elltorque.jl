
function get_dense_background(m::ModelSetup{T,D},les) where {T,D<:ModelDim}
    N,a,b,c,b0 = m.N,m.a,m.b,m.c,m.b0
    if D == Full
        LHS, RHS, vs1 = Mire.assemblemhd(N, a, b, c, [0,0,1/m.Le], b0)
        vs2 = vs1
    elseif D == Hybrid
        LHS, RHS, vs1, vs2 = Mire.assemblemhd_hybrid(N, N, a, b, c, [0,0,1/m.Le], b0)
    elseif D == QG
        LHS, RHS, vs1 = Mire.assemblemhd_qg(N, a, b, c, [0,0,1/m.Le], b0)
        vs2 = vs1
    elseif D == Geo
        LHS, RHS, vs1, vs2 = assemblemhd_geo(N, 2N+1, a, b, c, [0,0,1/m.Le], b0)
    else
        error("model must be one of: Full, Hybrid or QG!")
    end
    RHSt = copy(RHS)
    cmat = cacheint(2N,a,b,c)
    invL = inv(Matrix(LHS))
    ωs, us, bs, eks, ebs = [], [], [], [], []
     for le in les
        Mire.projectforce!(view(RHSt,1:length(vs2),1:length(vs2)),cmat,vs2,vs2,coriolis,[0,0,1/le])
        esol = eigen(Float64.(invL*RHSt))
        inds = 1e-2.<abs.(esol.values).<1e3
        u,b,ek,eb = Elltorque.get_ub(esol.vectors[:,inds],vs1,vs2,Float64.(cmat); ekin=false, getenergies=true)
        push!(ωs,esol.values[inds])
        push!(us,u)
        push!(bs,b)
        push!(eks,ek)
        push!(ebs,eb)
    end
    return ωs,us,bs,eks,ebs
end

function highres_avoidedcrossing(ω1,ω2,le0,le1,nle,m,LHS0,RHS0)
    les = range(le0,le1,length=nle)
    λs1,us1,lesout1 = tracking_lehnert(m, les, ω1, LHS0, RHS0,nev=1,verbose=false, maxdlefac=600);
    λs2,us2,lesout2 = tracking_lehnert(m, les, ω2, LHS0, RHS0,nev=1,verbose=false, maxdlefac=600);
    return λs1,us1,lesout1,λs2,us2,lesout2
end

a,b,c,Le = df64"1.25",df64"0.8",df64"1",df64"1e-5";

pAform = (y^0 + x)/(3*one(a))
b0Af = (a,b,c)-> b0_Aform(pAform,a,b,c)

if CALCULATE
    m = ModelSetup(a,b,c,Le,b0Af, "aform_ellipse2", 7, QG())
    LHS0,RHS0,vs_qg0 = assemblemhd_qg(m.N,m.a,m.b,m.c,[0,0,1/m.Le],m.b0)
    cmat = cacheint(m.N,m.a,m.b,m.c)
    esol=eigen(Float64.(inv(Matrix(LHS0))*RHS0))
    tms = esol.values[1e-1.<imag.(esol.values).<50]
    TARGETS = tms[sortperm(abs.(tms))]
    slowmodes = esol.values[abs.(esol.values).<1e-2]
    fastest_slowmode = slowmodes[sortperm(abs.(slowmodes))[end]]
    les=10.0.^range(-6,-2,length=200);

    λs1,us1,lesout1  = tracking_lehnert(m, les, TARGETS[1], LHS0, RHS0,nev=2,verbose=false)
    λs2,us2,lesout2 = tracking_lehnert(m, les, TARGETS[2], LHS0, RHS0,nev=2,verbose=false)
    λs3,us3,lesout3 = tracking_lehnert(m, les, TARGETS[3], LHS0, RHS0,nev=2,verbose=false)
    λs4,us4,lesout4 = tracking_lehnert(m, les, TARGETS[4], LHS0, RHS0,nev=2,verbose=false)
    λs5,us5,lesout5 = tracking_lehnert(m, les, fastest_slowmode, LHS0, RHS0,nev=2,verbose=false)

    v1,b1,ek1,eb1 = broadcast(us->Elltorque.get_ub1(us/mean(us),vs_qg0,vs_qg0,cmat,getenergies=true),us1) |> destruct;
    v2,b2,ek2,eb2 = broadcast(us->Elltorque.get_ub1(us/mean(us),vs_qg0,vs_qg0,cmat,getenergies=true),us2) |> destruct;
    v3,b3,ek3,eb3 = broadcast(us->Elltorque.get_ub1(us/mean(us),vs_qg0,vs_qg0,cmat,getenergies=true),us3) |> destruct;
    v4,b4,ek4,eb4 = broadcast(us->Elltorque.get_ub1(us/mean(us),vs_qg0,vs_qg0,cmat,getenergies=true),us4) |> destruct;
    v5,b5,ek5,eb5 = broadcast(us->Elltorque.get_ub1(us/mean(us),vs_qg0,vs_qg0,cmat,getenergies=true),us5) |> destruct;


    L1 = [Elltorque.angularmom(v,cmat,3) for v in v1];
    L2 = [Elltorque.angularmom(v,cmat,3) for v in v2];
    L3 = [Elltorque.angularmom(v,cmat,3) for v in v3];
    L4 = [Elltorque.angularmom(v,cmat,3) for v in v4];
    L5 = [Elltorque.angularmom(v,cmat,3) for v in v5];

    λs1_ac,us1_ac,lesout1_ac,λs2_ac,us2_ac,lesout2_ac = highres_avoidedcrossing(0.376im,0.355im,7.7e-4,8.3e-4,100,m,LHS0,RHS0);

    v1_ac,b1_ac,ek1_ac,eb1_ac = broadcast(us->Elltorque.get_ub1(us/mean(us),vs_qg0,vs_qg0,cmat,getenergies=true),us1_ac) |> destruct;
    v2_ac,b2_ac,ek2_ac,eb2_ac = broadcast(us->Elltorque.get_ub1(us/mean(us),vs_qg0,vs_qg0,cmat,getenergies=true),us2_ac) |> destruct;

    ωs,us,bs,eks,ebs = get_dense_background(m,les);

    usall = vcat(us...);
    bsall = vcat(bs...);
    angularmomall = broadcast(u->Elltorque.angularmom(u, cmat, 3), usall);
    lesall   = vcat([les[i]*ones(length(ωs[i])) for i=1:length(les)]...)
    evalsall = vcat([abs.(ωs[i]) for i=1:length(les)]...);
    ratiosall = vcat([abs.(eks[i])./abs.(ebs[i]) for i=1:length(les)]...);

    fname = joinpath(datapath, "tracking_lehnert_qg_mode1.jld2")
    JLD2.@save fname λs1 us1 lesout1 v1 L1 b1 ek1 eb1
    fname = joinpath(datapath, "tracking_lehnert_qg_mode2.jld2")
    JLD2.@save fname λs2 us2 lesout2 v2 L2 b2 ek2 eb2
    fname = joinpath(datapath, "tracking_lehnert_qg_mode3.jld2")
    JLD2.@save fname λs3 us3 lesout3 v3 L3 b3 ek3 eb3
    fname = joinpath(datapath, "tracking_lehnert_qg_mode4.jld2")
    JLD2.@save fname λs4 us4 lesout4 v4 L4 b4 ek4 eb4
    fname = joinpath(datapath, "tracking_lehnert_qg_mode5.jld2")
    JLD2.@save fname λs5 us5 lesout5 v5 L5 b5 ek5 eb5
    fname = joinpath(datapath, "tracking_lehnert_qg_dense_ratios.jld2")
    JLD2.@save fname lesall evalsall ratiosall usall bsall angularmomall
    fname = joinpath(datapath, "tracking_lehnert_qg_avoided_zoom.jld2")
    JLD2.@save fname  λs1_ac us1_ac lesout1_ac λs2_ac us2_ac lesout2_ac v1_ac b1_ac ek1_ac eb1_ac v2_ac b2_ac ek2_ac eb2_ac
else
    fname = joinpath(datapath, "tracking_lehnert_qg_mode1.jld2")
    JLD2.@load fname λs1 us1 lesout1 v1 L1 b1 ek1 eb1
    fname = joinpath(datapath, "tracking_lehnert_qg_mode2.jld2")
    JLD2.@load fname λs2 us2 lesout2 v2 L2 b2 ek2 eb2
    fname = joinpath(datapath, "tracking_lehnert_qg_mode3.jld2")
    JLD2.@load fname λs3 us3 lesout3 v3 L3 b3 ek3 eb3
    fname = joinpath(datapath, "tracking_lehnert_qg_mode4.jld2")
    JLD2.@load fname λs4 us4 lesout4 v4 L4 b4 ek4 eb4
    fname = joinpath(datapath, "tracking_lehnert_qg_mode5.jld2")
    JLD2.@load fname λs5 us5 lesout5 v5 L5 b5 ek5 eb5
    fname = joinpath(datapath, "tracking_lehnert_qg_dense_ratios.jld2")
    JLD2.@load fname lesall evalsall ratiosall usall bsall angularmomall
    fname = joinpath(datapath, "tracking_lehnert_qg_avoided_zoom.jld2")
    JLD2.@load fname  λs1_ac us1_ac lesout1_ac λs2_ac us2_ac lesout2_ac v1_ac b1_ac ek1_ac eb1_ac v2_ac b2_ac ek2_ac eb2_ac
end

PyPlot.rc("text",usetex=true)
# PyPlot.rc("figure", dpi=600)

f,ax=subplots(2,sharex=true,sharey=true,figsize=(3,3.5), dpi=600)

alphaval=0.35
lw=2.3

minratio = minimum(abs.(ek5./eb5)[7.85e-4.<lesout5.<8.05e-4])
maxratio = maximum(abs.(ek1./eb1)[7.7e-4.<lesout1.<8.2e-4])
minratio = minimum(ratiosall[1e-1 .< evalsall .<1e1])
maxratio = maximum(ratiosall[1e-1 .< evalsall .<1e1])

#axis 2

ax[2].loglog(range(1e-6,3e-3,length=20),abs(λs1[1])*ones(20),"-", color="C0", lw=lw,zorder=2)
ax[2].loglog(range(1e-6,4e-3,length=20),abs(λs2[1])*ones(20),"-.", color="C1", lw=lw,zorder=2)
ax[2].loglog(range(1e-6,4e-3,length=20),abs(λs3[1])*ones(20),"--", color="C2", lw=lw,zorder=2)
ax[2].loglog(range(1e-6,8e-3,length=20),abs(λs4[1])*ones(20),":", color="C3", lw=lw,zorder=2)
ax[2].loglog(range(1e-4,5e-3,length=10),
    abs.(λs5[lesout5.<1e-4][end])*range(1e-4,5e-3,length=10)/1e-4,
             ls=(0, (3, 1, 1, 1, 1, 1)),color="C4", lw=lw, zorder=1)

cm = :Greys

per = sortperm(ratiosall)
im2 = ax[2].scatter(lesall[per],evalsall[per], c = ratiosall[per], s=0.2,
    rasterized=true, norm=PyPlot.matplotlib.colors.LogNorm(vmin=minratio,vmax=maxratio), cmap=cm)

cb2 = colorbar(im2,ax=ax[2])
cb2.set_label(L"E_\mathrm{kin}/E_\mathrm{mag}")


#axis 1
per = sortperm(ratiosall)
im2 = ax[1].scatter(lesall[per],evalsall[per], c = ratiosall[per], s=0.2,
    rasterized=true, norm=PyPlot.matplotlib.colors.LogNorm(vmin=minratio,vmax=maxratio), cmap=cm)

cb2 = colorbar(im2,ax=ax[1])
cb2.set_label(L"E_\mathrm{kin}/E_\mathrm{mag}")



# k = 1:30:length(lesout1)
# k2 = 1:30:length(lesout5)


#inset:
axins = ax[1].inset_axes([0.1, 1.1, 0.8, 0.7])
im1=axins.scatter(lesout1_ac,abs.(λs1_ac), s=4,c=abs.(ek1_ac./eb1_ac),cmap=cm,norm=PyPlot.matplotlib.colors.LogNorm(vmin=minratio,vmax=maxratio))
axins.scatter(lesout2_ac,abs.(λs2_ac), s=4,c=abs.(ek2_ac./eb2_ac),cmap=cm,norm=PyPlot.matplotlib.colors.LogNorm(vmin=minratio,vmax=maxratio))

x1, x2, y1, y2 = 7.8e-4,8.15e-4, 0.37,0.382
axins.set_xlim([x1, x2])
axins.set_ylim([y1, y2])
axins.set_xticklabels("")
axins.set_yticklabels("")
ax[1].indicate_inset_zoom(axins)


#axis 1
# ax.set_xlim([1e-6,5e-2])
ax[1].set_ylim([1e-1,1e1])
# ""
# ax[1].set_xlabel(L"\mathrm{Le}")
ax[1].set_ylabel(L"|\omega|")

#axis 2
ax[2].set_xlim([1e-6,1e-2])
ax[2].set_ylim([1e-1,1e1])
# ""
ax[2].set_xlabel(L"\mathrm{Le}")
ax[2].set_ylabel(L"|\omega|")

if SAVEFIGS
    figname=joinpath(figpath,"avoidedcrossing_2.pdf")
    savefig(figname, bbox_inches="tight")
end





clf()






PyPlot.rc("text",usetex=true)
# PyPlot.rc("figure", dpi=600)
f,ax = subplots(2,sharex=true,figsize=(3,3.5), dpi=600)

alphaval = 0.15

ax[1].loglog(lesout1[lesout1.<7.5e-4],abs.(L1[lesout1.<7.5e-4]),"-", color="C0")
ax[1].loglog(lesout2[lesout2.<1.4e-3],abs.(L2[lesout2.<1.4e-3]),"-.",color="C1")
ax[1].loglog(lesout3[lesout3.<2.5e-3],abs.(L3[lesout3.<2.5e-3]),"--",color="C2")
ax[1].loglog(lesout4[lesout4.<4e-3],abs.(L4[lesout4.<4e-3]),":",color="C3")

ax[1].set_ylabel(L"|L_z|")

per = sortperm(abs.(evalsall))
clr=evalsall[per]
# clr=(0.4,0.4,0.4)

δω = 1/10
δL = 1/2
select=abs.(λs1[1])*(1-δω).<abs.(evalsall)[per].<abs.(λs1[1])*(1+δω)
select = select .& (lesall[per].>7.5e-4)
select = select .& (abs.(L1[1])*(1-δL).<abs.(angularmomall)[per].<abs.(L1[1])*(1+δL))
ax[1].scatter(lesall[per][select],abs.(angularmomall[per][select]), c = "C0", s=1,
    rasterized=true, norm=PyPlot.matplotlib.colors.LogNorm(vmin=1e-3,vmax=1e2))

select=abs.(λs2[1])*(1-δω).<abs.(evalsall)[per].<abs.(λs2[1])*(1+δω)
select = select .& (lesall[per].>1e-3)
select = select .& (abs.(L2[1])*(1-δL).<abs.(angularmomall)[per].<abs.(L2[1])*(1+δL))
ax[1].scatter(lesall[per][select],abs.(angularmomall[per][select]), c = "C1", s=1,
    rasterized=true, norm=PyPlot.matplotlib.colors.LogNorm(vmin=1e-3,vmax=1e2))

select=abs.(λs3[1])*(1-δω).<abs.(evalsall)[per].<abs.(λs3[1])*(1+δω)
select = select .& (lesall[per].>1e-3)
select = select .& (abs.(L3[1])*(1-δL).<abs.(angularmomall)[per].<abs.(L3[1])*(1+δL))
ax[1].scatter(lesall[per][select],abs.(angularmomall[per][select]), c = "C2", s=1,
    rasterized=true, norm=PyPlot.matplotlib.colors.LogNorm(vmin=1e-3,vmax=1e2))

select=abs.(λs4[1])*(1-δω).<abs.(evalsall)[per].<abs.(λs4[1])*(1+δω)
select = select .& (lesall[per].>1e-3)
select = select .& (abs.(L4[1])*(1-δL).<abs.(angularmomall)[per].<abs.(L4[1])*(1+δL))
ax[1].scatter(lesall[per][select],abs.(angularmomall[per][select]), c = "C3", s=1,
    rasterized=true, norm=PyPlot.matplotlib.colors.LogNorm(vmin=1e-3,vmax=1e2))

ax[1].set_ylim([1e-4,1])


ax[2].loglog(lesout1[lesout1.<7.5e-4],abs.(L1.*λs1)[lesout1.<7.5e-4],"-", color="C0")
ax[2].loglog(lesout2[lesout2.<1.4e-3],abs.(L2.*λs2)[lesout2.<1.4e-3],"-.",color="C1")
ax[2].loglog(lesout3[lesout3.<2.5e-3],abs.(L3.*λs3)[lesout3.<2.5e-3],"--",color="C2")
ax[2].loglog(lesout4[lesout4.<4e-3],abs.(L4.*λs4)[lesout4.<4e-3],":",color="C3")

δω = 1/10
select=abs.(λs1[1])*(1-δω).<abs.(evalsall)[per].<abs.(λs1[1])*(1+δω)
select = select .& (lesall[per].>7.5e-4)
select = select .& (abs.(L1[1])*(1-δL).<abs.(angularmomall)[per].<abs.(L1[1])*(1+δL))
ax[2].scatter(lesall[per][select],abs.(angularmomall.*evalsall)[per][select], c = "C0", s=1,
    rasterized=true, norm=PyPlot.matplotlib.colors.LogNorm(vmin=1e-3,vmax=1e2))

select=abs.(λs2[1])*(1-δω).<abs.(evalsall)[per].<abs.(λs2[1])*(1+δω)
select = select .& (lesall[per].>1e-3)
select = select .& (abs.(L2[1])*(1-δL).<abs.(angularmomall)[per].<abs.(L2[1])*(1+δL))
ax[2].scatter(lesall[per][select],abs.(angularmomall.*evalsall)[per][select], c = "C1", s=1,
    rasterized=true, norm=PyPlot.matplotlib.colors.LogNorm(vmin=1e-3,vmax=1e2))

select=abs.(λs3[1])*(1-δω).<abs.(evalsall)[per].<abs.(λs3[1])*(1+δω)
select = select .& (lesall[per].>1e-3)
select = select .& (abs.(L3[1])*(1-δL).<abs.(angularmomall)[per].<abs.(L3[1])*(1+δL))
ax[2].scatter(lesall[per][select],abs.(angularmomall.*evalsall)[per][select], c = "C2", s=1,
    rasterized=true, norm=PyPlot.matplotlib.colors.LogNorm(vmin=1e-3,vmax=1e2))

select=abs.(λs4[1])*(1-δω).<abs.(evalsall)[per].<abs.(λs4[1])*(1+δω)
select = select .& (lesall[per].>1e-3)
select = select .& (abs.(L4[1])*(1-δL).<abs.(angularmomall)[per].<abs.(L4[1])*(1+δL))
ax[2].scatter(lesall[per][select],abs.(angularmomall.*evalsall)[per][select], c = "C3", s=1,
    rasterized=true, norm=PyPlot.matplotlib.colors.LogNorm(vmin=1e-3,vmax=1e2))

ax[2].set_ylim([1e-4,2e-1])
ax[2].set_ylabel(L"|\omega L_z|")

ax[2].set_xlim([1e-6,1e-2])
ax[2].set_xlabel(L"\mathrm{Le}")

if SAVEFIGS
    figname=joinpath(figpath,"inertialtorque_z_le_qg_aform.pdf")
    savefig(figname, bbox_inches="tight")
end


if CALCULATE
    m = ModelSetup(a,b,c,Le,b0_2_8, "aform_ellipse2",9 , Hybrid())
    LHS0,RHS0,vs0,vs_qg0 = assemblemhd_hybrid(m.N,m.N,m.a,m.b,m.c,[0,0,1/m.Le],m.b0)
    cmat = cacheint(m.N,m.a,m.b,m.c)
    esol=eigen(Float64.(inv(Matrix(LHS0))*RHS0))
    tms = esol.values[1e-1.<imag.(esol.values).<50]
    TARGETS = tms[sortperm(abs.(tms))]
    slowmodes = esol.values[abs.(esol.values).<1e-2]
    fastest_slowmode = slowmodes[sortperm(abs.(slowmodes))[end]]

    les=10.0.^range(-6,-2,length=200);
    @show length(TARGETS)
    λs1,us1,lesout1  = tracking_lehnert(m, les, TARGETS[1], LHS0, RHS0,nev=2,verbose=false)
    λs2,us2,lesout2 = tracking_lehnert(m, les, TARGETS[2], LHS0, RHS0,nev=2,verbose=false)
    λs3,us3,lesout3 = tracking_lehnert(m, les, TARGETS[3], LHS0, RHS0,nev=2,verbose=false)
    λs4,us4,lesout4 = tracking_lehnert(m, les, TARGETS[4], LHS0, RHS0,nev=2,verbose=false)
    λs5,us5,lesout5 = tracking_lehnert(m, les, fastest_slowmode, LHS0, RHS0,nev=2,verbose=false);

    v1,b1,ek1,eb1 = broadcast(us->Elltorque.get_ub1(us/mean(us),vs0,vs_qg0,cmat,getenergies=true),us1) |> destruct;
    v2,b2,ek2,eb2 = broadcast(us->Elltorque.get_ub1(us/mean(us),vs0,vs_qg0,cmat,getenergies=true),us2) |> destruct;
    v3,b3,ek3,eb3 = broadcast(us->Elltorque.get_ub1(us/mean(us),vs0,vs_qg0,cmat,getenergies=true),us3) |> destruct;
    v4,b4,ek4,eb4 = broadcast(us->Elltorque.get_ub1(us/mean(us),vs0,vs_qg0,cmat,getenergies=true),us4) |> destruct;
    v5,b5,ek5,eb5 = broadcast(us->Elltorque.get_ub1(us/mean(us),vs0,vs_qg0,cmat,getenergies=true),us5) |> destruct;


    L1 = [Elltorque.angularmom(v,cmat,3) for v in v1];
    L2 = [Elltorque.angularmom(v,cmat,3) for v in v2];
    L3 = [Elltorque.angularmom(v,cmat,3) for v in v3];
    L4 = [Elltorque.angularmom(v,cmat,3) for v in v4];
    L5 = [Elltorque.angularmom(v,cmat,3) for v in v5];

    ωs,us,bs,eks,ebs = get_dense_background(m,les);

    usall = vcat(us...);
    bsall = vcat(bs...);
    angularmomall = broadcast(u->Elltorque.angularmom(u, cmat, 3), usall);
    lesall   = vcat([les[i]*ones(length(ωs[i])) for i=1:length(les)]...)
    evalsall = vcat([abs.(ωs[i]) for i=1:length(les)]...);
    ratiosall = vcat([abs.(eks[i])./abs.(ebs[i]) for i=1:length(les)]...);

    fname = joinpath(datapath, "tracking_lehnert_hyb_mode1.jld2")
    JLD2.@save fname λs1 us1 lesout1 v1 L1 b1 ek1 eb1
    fname = joinpath(datapath, "tracking_lehnert_hyb_mode2.jld2")
    JLD2.@save fname λs2 us2 lesout2 v2 L2 b2 ek2 eb2
    fname = joinpath(datapath, "tracking_lehnert_hyb_mode3.jld2")
    JLD2.@save fname λs3 us3 lesout3 v3 L3 b3 ek3 eb3
    fname = joinpath(datapath, "tracking_lehnert_hyb_mode4.jld2")
    JLD2.@save fname λs4 us4 lesout4 v4 L4 b4 ek4 eb4
    fname = joinpath(datapath, "tracking_lehnert_hyb_mode5.jld2")
    JLD2.@save fname λs5 us5 lesout5 v5 L5 b5 ek5 eb5
    fname = joinpath(datapath, "tracking_lehnert_hyb_dense.jld2")
    JLD2.@save fname lesall evalsall ratiosall usall bsall angularmomall
else
    fname = joinpath(datapath, "tracking_lehnert_hyb_mode1.jld2")
    JLD2.@load fname λs1 us1 lesout1 v1 L1 b1 ek1 eb1
    fname = joinpath(datapath, "tracking_lehnert_hyb_mode2.jld2")
    JLD2.@load fname λs2 us2 lesout2 v2 L2 b2 ek2 eb2
    fname = joinpath(datapath, "tracking_lehnert_hyb_mode3.jld2")
    JLD2.@load fname λs3 us3 lesout3 v3 L3 b3 ek3 eb3
    fname = joinpath(datapath, "tracking_lehnert_hyb_mode4.jld2")
    JLD2.@load fname λs4 us4 lesout4 v4 L4 b4 ek4 eb4
    fname = joinpath(datapath, "tracking_lehnert_hyb_mode5.jld2")
    JLD2.@load fname λs5 us5 lesout5 v5 L5 b5 ek5 eb5
    fname = joinpath(datapath, "tracking_lehnert_hyb_dense.jld2")
    JLD2.@load fname lesall evalsall ratiosall usall bsall angularmomall
end

PyPlot.rc("text",usetex=true)
PyPlot.rc("figure", dpi=600)
f,ax = subplots(2,sharex=true,figsize=(3,3.5))

alphaval = 0.15

ax[1].loglog(lesout1[lesout1.<7.5e-4],abs.(L1[lesout1.<7.5e-4]),"-", color="C0")
ax[1].loglog(lesout2[lesout2.<1.4e-3],abs.(L2[lesout2.<1.4e-3]),"-.",color="C1")
ax[1].loglog(lesout3[lesout3.<2.5e-3],abs.(L3[lesout3.<2.5e-3]),"--",color="C2")
ax[1].loglog(lesout4[lesout4.<4e-3],abs.(L4[lesout4.<4e-3]),":",color="C3")

ax[1].set_ylabel(L"|L_z|")

per = sortperm(abs.(evalsall))
clr=evalsall[per]
# clr=(0.4,0.4,0.4)

δω = 1/10
δL = 1/2
select=abs.(λs1[1])*(1-δω).<abs.(evalsall)[per].<abs.(λs1[1])*(1+δω)
select = select .& (lesall[per].>7.5e-4)
select = select .& (abs.(L1[1])*(1-δL).<abs.(angularmomall)[per].<abs.(L1[1])*(1+δL))
ax[1].scatter(lesall[per][select],abs.(angularmomall[per][select]), c = "C0", s=1,
    rasterized=true, norm=PyPlot.matplotlib.colors.LogNorm(vmin=1e-3,vmax=1e2))

select=abs.(λs2[1])*(1-δω).<abs.(evalsall)[per].<abs.(λs2[1])*(1+δω)
select = select .& (lesall[per].>1e-3)
select = select .& (abs.(L2[1])*(1-δL).<abs.(angularmomall)[per].<abs.(L2[1])*(1+δL))
ax[1].scatter(lesall[per][select],abs.(angularmomall[per][select]), c = "C1", s=1,
    rasterized=true, norm=PyPlot.matplotlib.colors.LogNorm(vmin=1e-3,vmax=1e2))

select=abs.(λs3[1])*(1-δω).<abs.(evalsall)[per].<abs.(λs3[1])*(1+δω)
select = select .& (lesall[per].>1e-3)
select = select .& (abs.(L3[1])*(1-δL).<abs.(angularmomall)[per].<abs.(L3[1])*(1+δL))
ax[1].scatter(lesall[per][select],abs.(angularmomall[per][select]), c = "C2", s=1,
    rasterized=true, norm=PyPlot.matplotlib.colors.LogNorm(vmin=1e-3,vmax=1e2))

select=abs.(λs4[1])*(1-δω).<abs.(evalsall)[per].<abs.(λs4[1])*(1+δω)
select = select .& (lesall[per].>1e-3)
select = select .& (abs.(L4[1])*(1-δL).<abs.(angularmomall)[per].<abs.(L4[1])*(1+δL))
ax[1].scatter(lesall[per][select],abs.(angularmomall[per][select]), c = "C3", s=1,
    rasterized=true, norm=PyPlot.matplotlib.colors.LogNorm(vmin=1e-3,vmax=1e2))
ax[1].set_ylim([1e-4,1])


ax[2].loglog(lesout1[lesout1.<7.5e-4],abs.(L1.*λs1)[lesout1.<7.5e-4],"-", color="C0")
ax[2].loglog(lesout2[lesout2.<1.4e-3],abs.(L2.*λs2)[lesout2.<1.4e-3],"-.",color="C1")
ax[2].loglog(lesout3[lesout3.<2.5e-3],abs.(L3.*λs3)[lesout3.<2.5e-3],"--",color="C2")
ax[2].loglog(lesout4[lesout4.<4e-3],abs.(L4.*λs4)[lesout4.<4e-3],":",color="C3")

δω = 1/10
select=abs.(λs1[1])*(1-δω).<abs.(evalsall)[per].<abs.(λs1[1])*(1+δω)
select = select .& (lesall[per].>7.5e-4)
select = select .& (abs.(L1[1])*(1-δL).<abs.(angularmomall)[per].<abs.(L1[1])*(1+δL))
ax[2].scatter(lesall[per][select],abs.(angularmomall.*evalsall)[per][select], c = "C0", s=1,
    rasterized=true, norm=PyPlot.matplotlib.colors.LogNorm(vmin=1e-3,vmax=1e2))

select=abs.(λs2[1])*(1-δω).<abs.(evalsall)[per].<abs.(λs2[1])*(1+δω)
select = select .& (lesall[per].>1e-3)
select = select .& (abs.(L2[1])*(1-δL).<abs.(angularmomall)[per].<abs.(L2[1])*(1+δL))
ax[2].scatter(lesall[per][select],abs.(angularmomall.*evalsall)[per][select], c = "C1", s=1,
    rasterized=true, norm=PyPlot.matplotlib.colors.LogNorm(vmin=1e-3,vmax=1e2))

select=abs.(λs3[1])*(1-δω).<abs.(evalsall)[per].<abs.(λs3[1])*(1+δω)
select = select .& (lesall[per].>1e-3)
select = select .& (abs.(L3[1])*(1-δL).<abs.(angularmomall)[per].<abs.(L3[1])*(1+δL))
ax[2].scatter(lesall[per][select],abs.(angularmomall.*evalsall)[per][select], c = "C2", s=1,
    rasterized=true, norm=PyPlot.matplotlib.colors.LogNorm(vmin=1e-3,vmax=1e2))

select=abs.(λs4[1])*(1-δω).<abs.(evalsall)[per].<abs.(λs4[1])*(1+δω)
select = select .& (lesall[per].>1e-3)
select = select .& (abs.(L4[1])*(1-δL).<abs.(angularmomall)[per].<abs.(L4[1])*(1+δL))
ax[2].scatter(lesall[per][select],abs.(angularmomall.*evalsall)[per][select], c = "C3", s=1,
    rasterized=true, norm=PyPlot.matplotlib.colors.LogNorm(vmin=1e-3,vmax=1e2))



ax[2].set_ylim([1e-4,4e-1])


ax[2].set_ylabel(L"|\omega L_z|")

ax[2].set_xlim([1e-6,1e-2])
ax[2].set_xlabel(L"\mathrm{Le}")

if SAVEFIGS
    figname=joinpath(figpath,"inertialtorque_z_le_hyb_b028.pdf")
    savefig(figname, bbox_inches="tight")
end
