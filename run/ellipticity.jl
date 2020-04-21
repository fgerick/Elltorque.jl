
a0,b0,c0,Le = 1.25,0.8,1.0,1e-7
eps0 = (a0^2-b0^2)/(a0^2+b0^2)

pAform = (y^0 + x)/3
b0Af = (a,b,c)-> b0_Aform(pAform,a,b,c)

m = ModelSetup(a0,b0,c0,Le,b0Af, "aform_ellipse", 7, QG());

LHS0,RHS0,vs0 = assemblemhd_qg(m.N,m.a,m.b,m.c,[0,0,1/Le],m.b0);

esol = eigen(Matrix(RHS0),Matrix(LHS0));

tms = esol.values[1e-2.<imag.(esol.values).<1e2];

TARGETS = tms[sortperm(abs.(tms))]

if CALCULATE
    λs1,us1,ϵsout1,vss1,cmats1 = tracking_ellipt_reverse(m, eps0 , 1e-3, TARGETS[1], LHS0, RHS0, b0Af; verbose=false);
    println("qg 1 done")
    λs2,us2,ϵsout2,vss2,cmats2 = tracking_ellipt_reverse(m, eps0 , 1e-3, TARGETS[2], LHS0, RHS0, b0Af; verbose=false);
    println("qg 2 done")
    λs3,us3,ϵsout3,vss3,cmats3 = tracking_ellipt_reverse(m, eps0 , 1e-3, TARGETS[3], LHS0, RHS0, b0Af; verbose=false);
    println("qg 3 done")
    λs4,us4,ϵsout4,vss4,cmats4 = tracking_ellipt_reverse(m, eps0 , 1e-3, TARGETS[4], LHS0, RHS0, b0Af; verbose=false);
    println("qg 4 done")

    v1,b1,ek1,eb1 = get_ub1.(us1./mean.(us1),vss1,vss1,cmats1,getenergies=true) |> destruct;
    println("qg ub 1 done")
    v2,b2,ek2,eb2 = get_ub1.(us2./mean.(us2),vss2,vss2,cmats2,getenergies=true) |> destruct;
    println("qg ub 2 done")
    v3,b3,ek3,eb3 = get_ub1.(us3./mean.(us3),vss3,vss3,cmats3,getenergies=true) |> destruct;
    println("qg ub 3 done")
    v4,b4,ek4,eb4 = get_ub1.(us4./mean.(us4),vss4,vss4,cmats4,getenergies=true) |> destruct;
    println("qg ub 4 done")

    L1=[angularmom(v1[i],cmats1[i]*pi,3) for i = 1:length(v1)];
    L2=[angularmom(v2[i],cmats2[i]*pi,3) for i = 1:length(v2)];
    L3=[angularmom(v3[i],cmats3[i]*pi,3) for i = 1:length(v3)];
    L4=[angularmom(v4[i],cmats4[i]*pi,3) for i = 1:length(v4)];
    println("qg amom done")

    fname = joinpath(datapath, "tracking_ellipse_qg_mode1.jld2")
    JLD2.@save fname λs1 us1 ϵsout1 vss1 cmats1 v1 b1 L1 ek1 eb1
    fname = joinpath(datapath, "tracking_ellipse_qg_mode2.jld2")
    JLD2.@save fname λs2 us2 ϵsout2 vss2 cmats2 v2 b2 L2 ek2 eb2
    fname = joinpath(datapath, "tracking_ellipse_qg_mode3.jld2")
    JLD2.@save fname λs3 us3 ϵsout3 vss3 cmats3 v3 b3 L3 ek3 eb3
    fname = joinpath(datapath, "tracking_ellipse_qg_mode4.jld2")
    JLD2.@save fname λs4 us4 ϵsout4 vss4 cmats4 v4 b4 L4 ek4 eb4
else
    fname = joinpath(datapath, "tracking_ellipse_qg_mode1.jld2")
    JLD2.@load fname λs1 us1 ϵsout1 vss1 cmats1 v1 b1 L1 ek1 eb1
    fname = joinpath(datapath, "tracking_ellipse_qg_mode2.jld2")
    JLD2.@load fname λs2 us2 ϵsout2 vss2 cmats2 v2 b2 L2 ek2 eb2
    fname = joinpath(datapath, "tracking_ellipse_qg_mode3.jld2")
    JLD2.@load fname λs3 us3 ϵsout3 vss3 cmats3 v3 b3 L3 ek3 eb3
    fname = joinpath(datapath, "tracking_ellipse_qg_mode4.jld2")
    JLD2.@load fname λs4 us4 ϵsout4 vss4 cmats4 v4 b4 L4 ek4 eb4
end

#Figure 8

PyPlot.rc("text",usetex=true)
f,ax = subplots(2,sharex=true,figsize=(3,3.5))

ax[1].loglog(ϵsout1,abs.(λs1),"-")
ax[1].loglog(ϵsout2,abs.(λs2),"-.")
ax[1].loglog(ϵsout3,abs.(λs3),"--")
ax[1].loglog(ϵsout4,abs.(λs4),":")
ax[1].annotate(L"\epsilon^{1/2}",(1e-3,5e-3))
ax[1].set_ylabel(L"|\omega|")
ax[1].set_ylim([1e-3,3])

ax[2].loglog(ϵsout1,abs.(λs1.-λs1[end]),"-")
ax[2].loglog(ϵsout2,abs.(λs2.-λs2[end]),"-.")
ax[2].loglog(ϵsout3,abs.(λs3.-λs3[end]),"--")
ax[2].loglog(ϵsout4,abs.(λs4.-λs4[end]),":")

ax[2].annotate(L"\epsilon^{1/2}",(1e-3,4e-2))
ax[2].annotate(L"\epsilon",(1e-3,1e-4))

ax[2].set_xlim([1e-4,eps0])
ax[2].set_ylim([1e-5,1])
ax[2].set_ylabel(L"|\omega(\epsilon)-\omega(\epsilon=0)|")
ax[2].set_xlabel(L"\epsilon")

if SAVEFIGS
    figname = joinpath(figpath,"frequency_ellipticity_qg.pdf")
    savefig(figname, bbox_inches="tight");
end


#Figure 10

PyPlot.rc("text",usetex=true)
f,ax = subplots(2,sharex=true,figsize=(3,3.5))


ax[1].loglog(ϵsout1,abs.(L1),"-")
ax[1].loglog(ϵsout2,abs.(L2),"-.")
ax[1].loglog(ϵsout3,abs.(L3),"--")
ax[1].loglog(ϵsout4,abs.(L4),":")

ax[1].annotate(L"\epsilon^{1/2}",(3e-5,7e-2))
ax[1].annotate(L"\epsilon",(1e-2,8e-6))

ax[1].set_ylabel(L"|L_z|")
ax[1].set_ylim([1e-7,5])

ax[2].loglog(ϵsout1,abs.(L1.*λs1),"-")
ax[2].loglog(ϵsout2,abs.(L2.*λs2),"-.")
ax[2].loglog(ϵsout3,abs.(L3.*λs3),"--")
ax[2].loglog(ϵsout4,abs.(L4.*λs4),":")

ax[2].annotate(L"\epsilon",(1e-2,8e-6))
ax[2].set_ylabel(L"|\omega L_z|")

ax[2].set_xlim([2e-5,eps0])
ax[2].set_ylim([1e-7,1])

ax[2].set_xlabel(L"\epsilon")

if SAVEFIGS
    figname = joinpath(figpath,"inertialtorque_ellipticity_qg_a.pdf")
    savefig(figname, bbox_inches="tight");
end



if CALCULATE

    m = ModelSetup(a0,b0,c0,Le,b0_2_8, "aform_b028", 9, Hybrid());
    LHS0,RHS0,vs0 = assemblemhd_hybrid(m.N,m.N,m.a,m.b,m.c,[0,0,1/Le],m.b0);

    esol = eigen(Matrix(RHS0),Matrix(LHS0));
    tms = esol.values[1e-2.<imag.(esol.values).<1e2];
    TARGETS = tms[sortperm(abs.(tms))]

    λs1,us1,ϵsout1,vss1,cmats1 = tracking_ellipt_reverse(m, eps0 , 1e-3, TARGETS[1], LHS0, RHS0, b0_2_8; zerothresh=1e-10, verbose=false);
    println("hyb 1 done")
    λs2,us2,ϵsout2,vss2,cmats2 = tracking_ellipt_reverse(m, eps0 , 1e-3, TARGETS[2], LHS0, RHS0, b0_2_8; verbose=false);
    println("hyb 2 done")
    λs3,us3,ϵsout3,vss3,cmats3 = tracking_ellipt_reverse(m, eps0 , 1e-3, TARGETS[3], LHS0, RHS0, b0_2_8; verbose=false);
    println("hyb 3 done")
    λs4,us4,ϵsout4,vss4,cmats4 = tracking_ellipt_reverse(m, eps0 , 1e-3, TARGETS[4], LHS0, RHS0, b0_2_8; verbose=false);
    println("hyb 4 done")
    v1,b1,ek1,eb1 = get_ub1.(us1./mean.(us1),getindex.(vss1,1),getindex.(vss1,2),cmats1,getenergies=true) |> destruct;
    println("hyb ub 1 done")
    v2,b2,ek2,eb2 = get_ub1.(us2./mean.(us2),getindex.(vss2,1),getindex.(vss2,2),cmats2,getenergies=true) |> destruct;
    println("hyb ub 2 done")
    v3,b3,ek3,eb3 = get_ub1.(us3./mean.(us3),getindex.(vss3,1),getindex.(vss3,2),cmats3,getenergies=true) |> destruct;
    println("hyb ub 3 done")
    v4,b4,ek4,eb4 = get_ub1.(us4./mean.(us4),getindex.(vss4,1),getindex.(vss4,2),cmats4,getenergies=true) |> destruct;
    println("hyb ub 4 done")
    L1=[angularmom(v1[i],cmats1[i]*pi,3) for i = 1:length(v1)];
    L2=[angularmom(v2[i],cmats2[i]*pi,3) for i = 1:length(v2)];
    L3=[angularmom(v3[i],cmats3[i]*pi,3) for i = 1:length(v3)];
    L4=[angularmom(v4[i],cmats4[i]*pi,3) for i = 1:length(v4)];
    println("hyb amom done")

    fname = joinpath(datapath, "tracking_ellipse_hyb_mode1.jld2")
    JLD2.@save fname λs1 us1 ϵsout1 vss1 cmats1 v1 b1 L1 ek1 eb1
    fname = joinpath(datapath, "tracking_ellipse_hyb_mode2.jld2")
    JLD2.@save fname λs2 us2 ϵsout2 vss2 cmats2 v2 b2 L2 ek2 eb2
    fname = joinpath(datapath, "tracking_ellipse_hyb_mode3.jld2")
    JLD2.@save fname λs3 us3 ϵsout3 vss3 cmats3 v3 b3 L3 ek3 eb3
    fname = joinpath(datapath, "tracking_ellipse_hyb_mode4.jld2")
    JLD2.@save fname λs4 us4 ϵsout4 vss4 cmats4 v4 b4 L4 ek4 eb4
else
    fname = joinpath(datapath, "tracking_ellipse_hyb_mode1.jld2")
    JLD2.@load fname λs1 us1 ϵsout1 vss1 cmats1 v1 b1 L1 ek1 eb1
    fname = joinpath(datapath, "tracking_ellipse_hyb_mode2.jld2")
    JLD2.@load fname λs2 us2 ϵsout2 vss2 cmats2 v2 b2 L2 ek2 eb2
    fname = joinpath(datapath, "tracking_ellipse_hyb_mode3.jld2")
    JLD2.@load fname λs3 us3 ϵsout3 vss3 cmats3 v3 b3 L3 ek3 eb3
    fname = joinpath(datapath, "tracking_ellipse_hyb_mode4.jld2")
    JLD2.@load fname λs4 us4 ϵsout4 vss4 cmats4 v4 b4 L4 ek4 eb4
end

#Figure A3

PyPlot.rc("text",usetex=true)
f,ax = subplots(2,sharex=true,figsize=(3,3.5))

ax[1].loglog(ϵsout1,abs.(λs1),"-")
ax[1].loglog(ϵsout2,abs.(λs2),"-.")
ax[1].loglog(ϵsout3,abs.(λs3),"--")
ax[1].loglog(ϵsout4,abs.(λs4),":")
ax[1].annotate(L"\epsilon^{1/2}",(1e-3,5e-3))
ax[1].set_ylabel(L"|\omega|")
ax[1].set_ylim([1e-3,5])

ax[2].loglog(ϵsout1,abs.(λs1.-λs1[end]),"-")
ax[2].loglog(ϵsout2,abs.(λs2.-λs2[end]),"-.")
ax[2].loglog(ϵsout3,abs.(λs3.-λs3[end]),"--")
ax[2].loglog(ϵsout4,abs.(λs4.-λs4[end]),":")

ax[2].annotate(L"\epsilon^{1/2}",(1e-3,4e-2))
ax[2].annotate(L"\epsilon",(1e-3,1e-4))

ax[2].set_xlim([1e-4,eps0])
ax[2].set_ylim([1e-5,1])
ax[2].set_ylabel(L"|\omega(\epsilon)-\omega(\epsilon=0)|")
ax[2].set_xlabel(L"\epsilon")

if SAVEFIGS
    figname = joinpath(figpath,"frequency_ellipticity_hyb.pdf")
    savefig(figname, bbox_inches="tight");
end


# Figure A4

PyPlot.rc("text",usetex=true)
f,ax = subplots(3,sharex=true,figsize=(3,5.25))


ax[1].loglog(ϵsout1,abs.(λs1),"-")
ax[1].loglog(ϵsout2,abs.(λs2),"-.")
ax[1].loglog(ϵsout3,abs.(λs3),"--")
ax[1].loglog(ϵsout4,abs.(λs4),":")
ax[1].annotate(L"\epsilon^{1/2}",(1e-3,5e-3))
ax[1].set_ylabel(L"|\omega|")
ax[1].set_ylim([1e-3,5])


ax[2].loglog(ϵsout1,abs.(L1),"-")
ax[2].loglog(ϵsout2,abs.(L2),"-.")
ax[2].loglog(ϵsout3,abs.(L3),"--")
ax[2].loglog(ϵsout4,abs.(L4),":")

ax[2].annotate(L"\epsilon^{1/2}",(3e-5,7e-2))
ax[2].annotate(L"\epsilon",(1e-2,4e-6))

ax[2].set_ylabel(L"|L_z|")
ax[2].set_ylim([1e-7,5])

ax[3].loglog(ϵsout1,abs.(L1.*λs1),"-")
ax[3].loglog(ϵsout2,abs.(L2.*λs2),"-.")
ax[3].loglog(ϵsout3,abs.(L3.*λs3),"--")
ax[3].loglog(ϵsout4,abs.(L4.*λs4),":")

ax[3].annotate(L"\epsilon",(1e-2,8e-6))
ax[3].set_ylabel(L"|\omega L_z|")

ax[3].set_xlim([2e-5,eps0])
ax[3].set_ylim([1e-7,1])

ax[3].set_xlabel(L"\epsilon")

if SAVEFIGS
    figname = joinpath(figpath,"inertialtorque_ellipticity_hyb.pdf")
    savefig(figname, bbox_inches="tight");
end


true
