# Tracking slow and fast Malkus modes as a function of ellipticity.

function track_malkus(SAVEDATA, datapath="")
    Le = 1e-8
    a,b,c=1.0,1.0,1.0

    N = 7
    N_3d = 7
    Ω = 1/Le*Mire.ez

    # targetsf_qg = [1.764705882350379e7im, 2.31e7im, 2.4999e7im]
    targetsf_qg = [1.764705882350379e7im, 2.307e7im, 2.4999e7im]
    #slow modes QG

    targetss_qg = [6.66667e-8im, 1.93333e-7im, 3.9e-7im];

    #fast and slow modes Mire (just use Zhang rossby modes + Malkus)
    targetsf_3d = [ωMalkusFast(m,Le,1) for m in 1:3]
    targetss_3d = [ωMalkusSlow(m,Le,1) for m in 1:3]

    #QG
    LHS_qg,RHS_qg,vs_qg = assemblemhd_qg(N,a,b,c,Ω,b0_1_3(a,b,c))
    #Hybrid
    LHS_hyb,RHS_hyb,vs,vs_qg = assemblemhd_hybrid(N,N,a,b,c,Ω,b0_1_3(a,b,c))
    #3D
    LHS_3d,RHS_3d,vs = assemblemhd(N,a,b,c,Ω,b0_1_3(a,b,c));

    m0_qg = ModelSetup(a,b,c,Le, b0_1_3,"malkus_qg",N, QG())
    m0_hyb = ModelSetup(a,b,c,Le, b0_1_3,"malkus_hyb",N, Hybrid())
    m0_3d = ModelSetup(a,b,c,Le, b0_1_3,"malkus_3d",N_3d, Full())



    # bs = range(1.,stop=0.7,length=30);
    ϵs = vcat([0.0],10.0.^range(-7,log10.(0.5),length=30));

    #fast modes
    imax=1000
    λ_qg_f2, u_qg_f2, ϵ_qg_f2 = tracking_ellipt(m0_qg,ϵs,targetsf_qg[2],LHS_qg,RHS_qg,b0_1_3; nev=2,maxiter=imax)
    println("qg2")
    flush(stdout)
    λ_qg_f3, u_qg_f3, ϵ_qg_f3 = tracking_ellipt(m0_qg,ϵs,targetsf_qg[3],LHS_qg,RHS_qg,b0_1_3; nev=2,maxiter=imax)
    println("qg3")
    flush(stdout)
    λ_hyb_f2, u_hyb_f2, ϵ_hyb_f2 = tracking_ellipt(m0_hyb,ϵs,targetsf_qg[2],LHS_hyb,RHS_hyb,b0_1_3; nev=2,maxiter=imax)
    println("hyb2")
    flush(stdout)
    λ_hyb_f3, u_hyb_f3, ϵ_hyb_f3 = tracking_ellipt(m0_hyb,ϵs,targetsf_qg[3],LHS_hyb,RHS_hyb,b0_1_3; nev=2,maxiter=imax)
    println("hyb3")
    flush(stdout)
    λ_3d_f2, u_3d_f2, ϵ_3d_f2 = tracking_ellipt(m0_3d,ϵs,targetsf_3d[2],LHS_3d,RHS_3d,b0_1_3; nev=2,maxiter=imax)
    println("3d2")
    flush(stdout)
    λ_3d_f3, u_3d_f3, ϵ_3d_f3 = tracking_ellipt(m0_3d,ϵs,targetsf_3d[3],LHS_3d,RHS_3d,b0_1_3; nev=2,maxiter=imax)
    println("3d3")
    flush(stdout)

    λ_qg_s2, u_qg_s2, ϵ_qg_s2 = tracking_ellipt(m0_qg,ϵs,targetss_qg[2],LHS_qg,RHS_qg,b0_1_3; nev=2,maxiter=imax)
    println("qg2s")
    flush(stdout)
    λ_qg_s3, u_qg_s3, ϵ_qg_s3 = tracking_ellipt(m0_qg,ϵs,targetss_qg[3],LHS_qg,RHS_qg,b0_1_3; nev=2,maxiter=imax)
    println("qg3s")
    flush(stdout)
    λ_hyb_s2, u_hyb_s2, ϵ_hyb_s2 = tracking_ellipt(m0_hyb,ϵs,targetss_qg[2],LHS_hyb,RHS_hyb,b0_1_3; nev=2,maxiter=imax)
    println("hyb2s")
    flush(stdout)
    λ_hyb_s3, u_hyb_s3, ϵ_hyb_s3 = tracking_ellipt(m0_hyb,ϵs,targetss_qg[3],LHS_hyb,RHS_hyb,b0_1_3; nev=2,maxiter=imax)
    println("hyb3s")
    flush(stdout)
    λ_3d_s2, u_3d_s2, ϵ_3d_s2 = tracking_ellipt(m0_3d,ϵs,1.9246950765951556e-7im,LHS_3d,RHS_3d,b0_1_3; nev=2,maxiter=imax)
    println("3d2s")
    flush(stdout)
    λ_3d_s3, u_3d_s3, ϵ_3d_s3 = tracking_ellipt(m0_3d,ϵs,targetss_3d[3],LHS_3d,RHS_3d,b0_1_3; nev=2,maxiter=imax)
    println("3d3s")
    flush(stdout)

    if SAVEDATA
        fname = joinpath(datapath,"malkus_ellipt_fast_qg.jld")
        JLD2.@save fname λ_qg_f2 u_qg_f2 ϵ_qg_f2 λ_qg_f3 u_qg_f3 ϵ_qg_f3

        fname = joinpath(datapath,"malkus_ellipt_fast_hyb.jld")
        JLD2.@save fname λ_hyb_f2 u_hyb_f2 ϵ_hyb_f2 λ_hyb_f3 u_hyb_f3 ϵ_hyb_f3

        fname = joinpath(datapath,"malkus_ellipt_fast_3d.jld")
        JLD2.@save fname λ_3d_f2 u_3d_f2 ϵ_3d_f2 λ_3d_f3 u_3d_f3 ϵ_3d_f3

        fname = joinpath(datapath,"malkus_ellipt_slow_qg.jld")
        JLD2.@save fname λ_qg_s2 u_qg_s2 ϵ_qg_s2 λ_qg_s3 u_qg_s3 ϵ_qg_s3

        fname = joinpath(datapath,"malkus_ellipt_slow_hyb.jld")
        JLD2.@save fname λ_hyb_s2 u_hyb_s2 ϵ_hyb_s2 λ_hyb_s3 u_hyb_s3 ϵ_hyb_s3

        fname = joinpath(datapath,"malkus_ellipt_slow_3d.jld")
        JLD2.@save fname λ_3d_s2 u_3d_s2 ϵ_3d_s2 λ_3d_s3 u_3d_s3 ϵ_3d_s3

    end
    return true
end


function plotmalkusbenchmark(datapath,figpath,SAVEFIGS)
    fname = joinpath(datapath,"malkus_ellipt_fast_qg.jld")
    JLD2.@load fname λ_qg_f2 u_qg_f2 ϵ_qg_f2 λ_qg_f3 u_qg_f3 ϵ_qg_f3

    fname = joinpath(datapath,"malkus_ellipt_fast_hyb.jld")
    JLD2.@load fname λ_hyb_f2 u_hyb_f2 ϵ_hyb_f2 λ_hyb_f3 u_hyb_f3 ϵ_hyb_f3

    fname = joinpath(datapath,"malkus_ellipt_fast_3d.jld")
    JLD2.@load fname λ_3d_f2 u_3d_f2 ϵ_3d_f2 λ_3d_f3 u_3d_f3 ϵ_3d_f3

    fname = joinpath(datapath,"malkus_ellipt_slow_qg.jld")
    JLD2.@load fname λ_qg_s2 u_qg_s2 ϵ_qg_s2 λ_qg_s3 u_qg_s3 ϵ_qg_s3

    fname = joinpath(datapath,"malkus_ellipt_slow_hyb.jld")
    JLD2.@load fname λ_hyb_s2 u_hyb_s2 ϵ_hyb_s2 λ_hyb_s3 u_hyb_s3 ϵ_hyb_s3

    fname = joinpath(datapath,"malkus_ellipt_slow_3d.jld")
    JLD2.@load fname λ_3d_s2 u_3d_s2 ϵ_3d_s2 λ_3d_s3 u_3d_s3 ϵ_3d_s3

    ell(a,b) = (a^2-b^2)/(a^2+b^2)

 # f,ax = subplots(nrows=2,ncols=1,sharex=true,figsize=(4,6));
    f,ax = subplots(2,sharex=true, figsize=(3,5))
    #fast modes subplot
    ax[1].plot(ϵ_qg_f2,abs.(λ_qg_f2), "k-", label="QG")
    ax[1].plot(ϵ_hyb_f2,abs.(λ_hyb_f2), "k.", label="Hybrid",markersize=5)
    ax[1].plot(ϵ_3d_f2,abs.(λ_3d_f2), "k--", label="3D")

    ax[1].set_ylabel(L"|\omega_{f,m=2}|")
    ax[1].set_xscale("log")
    ax[1].set_ylim([2.05e7,2.35e7])


    ax2=ax[1].twinx()
    ax2.plot(ϵ_qg_f3,abs.(λ_qg_f3), "C0-", label="QG")
    ax2.plot(ϵ_hyb_f3,abs.(λ_hyb_f3), "C0.", label="Hybrid",markersize=5)
    ax2.plot(ϵ_3d_f3,abs.(λ_3d_f3), "C0--", label="3D")

    ax2.set_ylim([2.35e7,2.65e7])
    ax2.set_ylabel(L"|\omega_{f,m=3}|",color="C0")
    ax2.tick_params("y",colors="C0")
    ax2.set_xscale("log")

    # ax[1].set_xlabel(L"\epsilon")
    # ax[1].set_xlim([1e-2,0.45])

    # if SAVEFIGS
    #     savefig(joinpath(figpath,"fastmodes_m23.pdf"),bbox_inches="tight")
    # end

    # #slow modes subplot

     # f,ax = subplots(1,figsize=(3,2))


    ax[2].plot(ϵ_qg_s2,abs.(λ_qg_s2), "k-", label="QG")
    ax[2].plot(ϵ_hyb_s2,abs.(λ_hyb_s2), "k.", label="Hybrid", markersize=5)
    ax[2].plot(ϵ_3d_s2,abs.(λ_3d_s2), "k--", label="3D")

    ax[2].set_ylabel(L"|\omega_{s,m=2}|")
    ax[2].set_ylim([1.9e-7,2.3e-7])

    ax3=ax[2].twinx()
    ax3.plot(ϵ_qg_s3,abs.(λ_qg_s3), "C0-", label="QG")
    ax3.plot(ϵ_hyb_s3,abs.(λ_hyb_s3), "C0.", label="Hybrid", markersize=5)
    ax3.plot(ϵ_3d_s3,abs.(λ_3d_s3), "C0--", label="3D")

    ax3.set_xscale("log")
    ax3.set_ylim([3.6e-7,4.6e-7])
    ax3.set_ylabel(L"|\omega_{s,m=3}|",color="C0")
    ax3.tick_params("y",colors="C0")


    # # ax[2].set_ylabel("|ω|")
    # # ax[2].set_ylim([3e-8,4.6e-7])
    # ax[2].legend(bbox_to_anchor=(0.4,1.)) #,bbox_to_anchor=(0.0,.91))
    ax[2].set_xlabel(L"\epsilon")
    ax[2].set_xlim([1e-2,0.45])

    if SAVEFIGS
        savefig(joinpath(figpath,"malkusmodes_m23.pdf"),bbox_inches="tight")
    end
end
