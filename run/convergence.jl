# calculate the convergence of the TM in the QG model (Figure 4)

if CALCULATE

    function torsion_2D(n,a,b,c,Ω,asymmfield,thresh_l,thresh_u,cmat; verbose=true)
        LHS,RHS,vs_qg = Mire.assemblemhd_qg(n,a,b,c,Ω,b0,cmat=cmat)
        evals = eigvals(Matrix(inv(Matrix(LHS))*RHS));
        torsionals = eachindex(evals)[(thresh_l .< abs.(imag.(evals)) .< thresh_u)]
        if verbose
            println("N=$n :")
            println(imag.(evals[torsionals]))
        end
        return evals[torsionals]
    end

    function b0_Aform(p_xy,a,b,c)
        ts = terms(p_xy)
        out = zero(Mire.qg_vel(0,0,a,b,c))
        for t in ts
            out .+= coefficient(t).*Mire.qg_vel(exponents(t)...,a,b,c)
        end
        return out
    end

    a,b,c = df64"1.25",df64"0.8", df64"1."
    Le = df64"1e-8"
    Ω = [0,0,1/Le]
    b0 = b0_Aform((y^0+x)/df64"3",a,b,c)
    ns = 1:2:25

    cmat = cacheint(27,a,b,c)
    tmodesqg_df64 = [torsion_2D(N,a,b,c,Ω,b0,1e-2,100.,cmat,verbose=VERBOSE) for N in ns]
    fname = joinpath(datapath,"qg_torsional_conversion.jld2")
    JLD2.@save fname tmodesqg_df64 ns
else
    fname = joinpath(datapath,"qg_torsional_conversion.jld2")
    JLD2.@load fname tmodesqg_df64 ns
end

#plotting:

figure(figsize=(4,3))
PyPlot.rc("text",usetex=true)
tmodspos = [t[imag.(t).>0] for t in tmodesqg_df64]
tsortper = sortperm.(imag.(tmodspos))
tsorted = [t[per] for (t,per) in zip(tmodspos,tsortper)]

for itorsion = 1:length(ns)-1
    if itorsion==1 #U_3 mode in gray
        C="gray"
        mark = "^"
    else
        C="k"
        mark="."
    end
    plot(ns[itorsion:end],getindex.(imag.(tsorted)[itorsion:end],itorsion),color=C,linestyle="-",marker=mark)
end
xlim([1,25])
ylim([0,6])
xticks(ns)
xlabel(L"N")
ylabel(L"|\omega|")

if SAVEFIGS
    figname = joinpath(figpath,"qg_convergence.pdf")
    savefig(figname,bbox_inches="tight");
end


true
