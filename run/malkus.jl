if CALCULATE
    track_malkus(true, datapath)
end

function plotit(SAVEFIGS::Bool)



      # diff frequencies:
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


      clf()
       f,ax = subplots(3,sharex=true, figsize=(3,7))
       #fast modes subplot
       ax[1].plot(ϵ_qg_f2,abs.(λ_qg_f2), "k-", label="QG")
       ax[1].plot(ϵ_hyb_f2,abs.(λ_hyb_f2), "k.", label="Hybrid",markersize=5)
       ax[1].plot(ϵ_3d_f2,abs.(λ_3d_f2), "k--", label="3D")

       ax[1].set_ylabel(L"|\omega_{f,m=2}|",color="k")
       ax[1].tick_params("y",colors="k")
       ax[1].set_xscale("log")
       ax[1].set_ylim([2.05e7,2.35e7])


       ax2=ax[1].twinx()
       c2="C0"
       ax2.plot(ϵ_qg_f3,abs.(λ_qg_f3), c2*"-", label="QG")
       ax2.plot(ϵ_hyb_f3,abs.(λ_hyb_f3), c2*".", label="Hybrid",markersize=5)
       ax2.plot(ϵ_3d_f3,abs.(λ_3d_f3), c2*"--", label="3D")
       ax2.set_ylim([2.35e7,2.65e7])
       ax2.set_ylabel(L"|\omega_{f,m=3}|",color=c2)
       ax2.tick_params("y",colors=c2)
       ax2.set_xscale("log")



       ax[2].plot(ϵ_qg_s2,abs.(λ_qg_s2), color="Gray", linestyle="-", label="QG")
       ax[2].plot(ϵ_hyb_s2,abs.(λ_hyb_s2), color="Gray", linestyle="",marker=".", label="Hybrid", markersize=5)
       ax[2].plot(ϵ_3d_s2,abs.(λ_3d_s2), color="Gray", linestyle="--", label="3D")

       ax[2].set_ylabel(L"|\omega_{s,m=2}|",color="Gray")
       ax[2].set_ylim([1.9e-7,2.3e-7])
       ax[2].tick_params("y",colors="Gray")

       ax3=ax[2].twinx()
       c4="C1"
       ax3.plot(ϵ_qg_s3,abs.(λ_qg_s3), c4*"-", label="QG")
       ax3.plot(ϵ_hyb_s3,abs.(λ_hyb_s3), c4*".", label="Hybrid", markersize=5)
       ax3.plot(ϵ_3d_s3,abs.(λ_3d_s3), c4*"--", label="3D")

       ax3.set_xscale("log")
       ax3.set_ylim([3.6e-7,4.6e-7])
       ax3.set_ylabel(L"|\omega_{s,m=3}|",color=c4)
       ax3.tick_params("y",colors=c4)

     ax[3].loglog(ϵ_qg_f2,abs.(λ_qg_f2[1].-λ_qg_f2)./abs.(λ_qg_f2[1]),"k-")
     ax[3].loglog(ϵ_hyb_f2,abs.(λ_hyb_f2[1].-λ_hyb_f2)./abs.(λ_hyb_f2[1]),"k.", markersize=5)
     ax[3].loglog(ϵ_3d_f2,abs.(λ_3d_f2[1].-λ_3d_f2)./abs.(λ_qg_f2[1]),"k--")

       ax[3].loglog(ϵ_qg_s3,abs.(λ_qg_s3[1].-λ_qg_s3)./abs.(λ_qg_s3[1]),"C0-")
       ax[3].loglog(ϵ_hyb_s3,abs.(λ_hyb_s3[1].-λ_hyb_s3)./abs.(λ_hyb_s3[1]),"C0.", markersize=5)
       ax[3].loglog(ϵ_3d_s3,abs.(λ_3d_s3[1].-λ_3d_s3)./abs.(λ_qg_s3[1]),"C0--")

       ax[3].loglog(ϵ_qg_s2,abs.(λ_qg_s2[1].-λ_qg_s2)./abs.(λ_qg_s2[1]), color="Gray", linestyle="-", label="QG")
       ax[3].loglog(ϵ_hyb_s2,abs.(λ_hyb_s2[1].-λ_hyb_s2)./abs.(λ_hyb_s2[1]), color="Gray", linestyle="",marker=".", label="Hybrid", markersize=5)
       ax[3].loglog(ϵ_3d_s2,abs.(λ_3d_s2[1].-λ_3d_s2)./abs.(λ_qg_s2[1]), color="Gray", linestyle="--", label="3D")

       ax[3].loglog(ϵ_qg_f3,abs.(λ_qg_f3[1].-λ_qg_f3)./abs.(λ_qg_f3[1]),"C1-")
       ax[3].loglog(ϵ_hyb_f3,abs.(λ_hyb_f3[1].-λ_hyb_f3)./abs.(λ_hyb_f3[1]),"C1.", markersize=5)
       ax[3].loglog(ϵ_3d_f3,abs.(λ_3d_f3[1].-λ_3d_f3)./abs.(λ_qg_f3[1]),"C1--")

       ax[3].annotate(L"\epsilon^{2}",(1e-1,5e-4))

       ax[3].set_ylabel(L"|\omega(\epsilon=0)-\omega|/|\omega(\epsilon=0)|")
       ax[3].set_xlabel(L"\epsilon")
       ax[3].set_xlim([1e-2,0.45])
       ax[3].set_ylim([1e-4,0.45^2])
       if SAVEFIGS
           savefig(joinpath(figpath,"malkusmodes_m23_v2.pdf"),bbox_inches="tight")
       end

       gcf()
    end

plotit(SAVEFIGS)
