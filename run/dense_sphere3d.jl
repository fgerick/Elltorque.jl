#function to get the velocities, magnetic field vectors and their kinetic and magnetic energy
#from the eigensolution.
function get_ub_energies(evecs,vs,vs_qg,cmat; ekin=false)
    nev=size(evecs,2)
    nqg=length(vs_qg)
    us,bs,ek,eb=[],[],[],[]
    for i=1:nev
        u = Mire.eigenvel(vs_qg,evecs[1:nqg,i])#/mean(evecs[:,i]))
        b = Mire.eigenvel(vs,evecs[nqg+1:end,i])#/mean(evecs[:,i]))
        energyu = Mire.inner_product(cmat,u,u)
        energyb = Mire.inner_product(cmat,b,b)
        push!(us,u)
        push!(bs,b)
        push!(ek,energyu)
        push!(eb,energyb)
    end
    return us,bs,ek,eb
end

#assemble matrices, calculate eigensolutions and their energies.
function calculate_energyratio(m::ModelSetup{T,D},datapath="",SAVEDATA=false,dtypename="f64") where {T,D <: ModelDim}
    a, b, c, Le, b0, N = m.a, m.b, m.c, m.Le, m.b0, m.N
    Ω = [0,0,1/Le]

    # different assembly function for different models, Full=3D, hybrid QG and 3D, and fully QG.
    if D==Full
        LHS, RHS, vs = Mire.assemblemhd(N, a, b, c, Ω, b0)
    elseif D==Hybrid
        LHS, RHS, vs, vs_qg = Mire.assemblemhd_hybrid(N, N, a, b, c, Ω, b0)
    elseif D==QG
        LHS, RHS, vs_qg = Mire.assemblemhd_qg(N, a, b, c, Ω, b0)
    else
        error("model must be one of: Full, Hybrid or QG!")
    end

    #cache the monomial integration
    cmat = cacheint(N, a, b, c)

    #solve the eigen problem
    if T <: LinearAlgebra.BlasFloat
        A,B = Matrix(RHS), Matrix(LHS)
        S = eigen(A, B)
    else
        A, B = complex.(Matrix(RHS)), complex.(Matrix(LHS))
        S = GenericSchur.schur(A, B)
    end

    #get eigenfrequencies and vectors
    ω = S.values
    evecs = eigvecs(S)

    if D == Full
        us, bs, ek, eb = get_ub_energies(evecs, vs, vs, cmat)
    elseif D == Hybrid
        us, bs, ek, eb = get_ub_energies(evecs, vs, vs_qg, cmat)
    elseif D == QG
        us, bs, ek, eb = get_ub_energies(evecs, vs_qg, vs_qg, cmat)
    end

    return ω, ek, eb
end


m = ModelSetup(1.0,1.0,1.0,1e-5, (a,b,c)->b0_1_3(a,b,c)+b0_1_1(a,b,c)/10,"sphere_fb1",5, Full());


a, b, c, Le, b0, N = m.a, m.b, m.c, m.Le, m.b0, m.N
Ω = [0,0,1/Le];

if CALCULATE
    @time LHS, RHS, vs = Mire.assemblemhd(N, a, b, c, Ω, b0);
    @time ω,ek,eb = calculate_energyratio(m);
    ωs,eks,ebs = [],[],[]
    les = 10.0.^range(-5,-1,length=50);
    for le in les
        println(le)
        flush(stdout)
        m = ModelSetup(1.0,1.0,1.0,le, (a,b,c)->b0_1_3(a,b,c)+b0_1_1(a,b,c)/10,"sphere_fb1",5, Full())
        ω,ek,eb = calculate_energyratio(m)
        push!(ωs,ω)
        push!(eks,ek)
        push!(ebs,eb)
    end
    fname = joinpath(datapath,"densespectrumsphere.jld2")
    JLD2.@save fname ωs eks ebs
else
    fname = joinpath(datapath,"densespectrumsphere.jld2")
    JLD2.@load fname ωs eks ebs
    les = 10.0.^range(-5,-1,length=50);
    end;

les_all = vcat([ones(length(eks[1]))*le for le in les]...)
evals_all = vcat(ωs...)
eks_all = vcat(eks...)
ebs_all = vcat(ebs...)
ratios_all = abs.(eks_all)./abs.(ebs_all);

perm = sortperm(ratios_all);


#Plot Figure 1 in paper:

PyPlot.rc("text",usetex=true)
f,ax=subplots(1,figsize=(4,3), dpi=600)


im2 = ax.scatter(les_all[perm],(abs.(imag.(evals_all)))[perm], c=ratios_all[perm],
                    norm=PyPlot.matplotlib.colors.LogNorm(vmin=1e-9,vmax=1e9),cmap=:cividis,s=8, rasterized=true);
cb = f.colorbar(im2,ax=ax)
cb.set_label(L"E_\mathrm{kin}/E_\mathrm{mag}")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim([les[1],les[end]])
ax.set_ylim([1e-6,1e6]);
ax.set_ylabel(L"|\omega|\, \left[ T_\Omega^{-1}\right]")
ax.set_ylabel(L"|\omega|")
ax.set_xlabel(L"\mathrm{Le}");


ax.annotate(L"\mathrm{slow} \propto \mathrm{Le}",(1e-3,5e-6),color="k")
ax.annotate(L"\mathrm{fast} \propto \mathrm{Le}^{-1}",(1e-3,1e4),color="k")
ax.annotate(L"\mathrm{TM}\propto \mathrm{const.}",(1.5e-5,0.85),color="k")

if SAVEFIGS
    figname=joinpath(figpath,"dense_sphere.pdf")
    savefig(figname,bbox_inches="tight")
end
