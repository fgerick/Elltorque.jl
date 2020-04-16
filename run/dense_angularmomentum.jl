

datapathn11 = "../data/runlehnert_b0An11"

if CALCULATE
    #This calculation takes a lot of time if run in serial. We recommend using
    #as many cores as available on one node, as the problem is embarissingly
    #parallel.
    np = 4 #number of cores
    Base.run(`julia dense_n11.jl $np $datapathn11`)
end

function loadall(m0::ModelSetup{T,D},Les,datapath,dtypename="df64") where {T,D<:ModelDim}

    ωso = []
    Lso = []
   for i=1:length(Les)
        m=ModelSetup{T,D}(m0.a,m0.b,m0.c,Les[i],m0.b0,"le_$i",m0.N)
        fname = joinpath(datapath,string(D)*"_$(m.name)_"*dtypename*"_N$(m.N).jld")
        fnameL = joinpath(datapath,"L_"*string(D)*"_$(m.name)_"*dtypename*"_N$(m.N).jld")
        if D == Full
            JLD2.@load fname ω
            JLD2.@load fnameL Ls
        elseif D == Hybrid
            JLD2.@load fname ω
            JLD2.@load fnameL Ls
        elseif D == QG
            JLD2.@load fname ω
            JLD2.@load fnameL Ls
        end
        push!(ωso,ω)
        push!(Lso,Ls)

    end

    return ωso, Lso
end



Les=df64"10.0".^range(-7,0,length=200);

m= ModelSetup(a,b,c,Double64(Les[1]),b0Af, "aform_ellipse_n11", 11, QG())

ωs,Ls = loadall(m,Les,datapathn11,"df64")
#all values in one vector:
evs_all = vcat(ωs...)
Ls_all = vcat(Ls...);
les_all = vcat([ones(length(ωs[1]))*le for le in Les]...);



#plotting:

PyPlot.rc("text",usetex=true)
figure(figsize=(4,3),dpi=600)
per = sortperm(abs.(Ls_all),rev=false)
per = per[imag.(evs_all[per]).>0]
scatter(les_all[per],abs.(evs_all)[per],c = abs.(evs_all[per].*Ls_all[per]), s = 1, rasterized = true,
    norm = PyPlot.matplotlib.colors.LogNorm(vmin = 1e-10,vmax = maximum(abs.(evs_all.*Ls_all))))
yscale("log")
xscale("log")
xlim([1e-6,5e-2])
ylim([1e-4,1e4])
cb = colorbar()
cb.set_label(L"|\omega L_z|")
xlabel(L"\mathrm{Le}")
ylabel(L"|\omega|")

if SAVEFIGS
    fname = joinpath(figpath,"dense_angularmom_n11.pdf")
    savefig(fname,bbox_inches="tight")
end
