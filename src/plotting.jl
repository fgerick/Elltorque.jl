function plottingtorquebalance(m::ModelSetup{T,D}, SAVEFIG=false) where {T<:Real,D<:ModelDim}
    PyPlot.rc("text",usetex=true)
    f,ax=subplots(3,sharex=true,figsize=(4,6));
    cstring = [L"\mathbf{\Gamma}_x",L"\mathbf{\Gamma}_y",L"\mathbf{\Gamma}_z"]
    for i=1:3
        figure()
        ax[i].plot(abs.(imag.(ω)), abs.(Γem[i])   , linewidth=0, color="C0", marker="s", markersize=10, label = L"|\mathbf{\Gamma}_\mathrm{b}|") #"|∫((b⋅∇)b₀+(b₀⋅∇)b)×r dV|")
        ax[i].plot(abs.(imag.(ω)), abs.(Γpmag[i]) , linewidth=0, color="C7", marker="d", markersize=9,  label = L"|\mathbf{\Gamma}_\mathrm{pm}|") #"|∫b⋅b₀(n×r) dS|")
    #     ax[i].plot(abs.(imag.(ω)), abs.(Γpageo[i]), linewidth=0, color="C8", marker="*", markersize=8,  label = "|∫(∇pₐ×r) dV|")
        ax[i].plot(abs.(imag.(ω)), abs.(Γp[i]) ,    linewidth=0, color="C2", marker=".", markersize=8,  label = L"|\mathbf{\Gamma}_\mathrm{p}|") #"|∫(∇p×r) dV|")
        ax[i].plot(abs.(imag.(ω)), abs.(Lω[i]) ,    linewidth=0, color="C1", marker="<", markersize=7,  label = L"|\mathbf{L}\omega|") #"|ω∫u×r dV|")
        ax[i].plot(abs.(imag.(ω)), abs.(Γptot[i]) , linewidth=0, color="C3", marker="+", markersize=6,  label = L"|\mathbf{\Gamma}_\mathrm{p}+\Gamma_\mathrm{pm}|")#"|∫∇(p+b⋅b₀)×r dV|")
        ax[i].plot(abs.(imag.(ω)), abs.(Γcor[i]),   linewidth=0, color="k",  marker="4", markersize=5,  label = L"|\mathbf{\Gamma}_{c}|") #"-2∫(Ω×u)×r dV")
    #     title(cstring[i]*" - coordinate")
        ax[i].set_ylabel(cstring[i])
        ax[i].set_yscale("symlog", linthreshy=1e-13)
        ax[i].set_yticks(vcat([0],10.0.^(-13:5:4)))
        ax[i].set_ylim([0,1e4])

    end
    ax[1].legend(framealpha=0.65, ncol=3, columnspacing=0.01)
    ax[3].set_xlabel(L"|\mathrm{Im}(\omega)|")
    ax[3].set_xscale("symlog", linthreshx=1e-6)
    ax[3].set_xticks(vcat([0],10.0.^(-6:3:6)))

    ax[3].set_xlim([0,1e6])
    f.savefig("figs/torquebalance_hybrid_Aform_n5.pdf", bbox_inches="tight")
