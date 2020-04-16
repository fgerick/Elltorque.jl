
function grid(x,y)
    X = [i for j in y, i in x]
    Y = [j for j in y, i in x]
    return X,Y
end

function plot_velocity_equator_uphi(a,b,v1; cbar=false, vmax=0, ngrid=50, kwargs...)
    X, Y = grid(range(-a,stop=a,length=ngrid),range(-b,stop=b,length=ngrid))
    ux =  real.([v1[1](Mire.x=>xt,Mire.y=>yt,Mire.z=>0) for (xt,yt) in zip(X,Y)])
    uy =  real.([v1[2](Mire.x=>xt,Mire.y=>yt,Mire.z=>0) for (xt,yt) in zip(X,Y)])
    radius = .√(X.^2/a^2+Y.^2/b^2)
    u = .√(ux.^2+uy.^2)
    phi = atan.(Y/b,X/a);
    uphi = @. (-a*ux*sin(phi) + b*uy*cos(phi))
    outsideellipse=(X.^2/a^2+Y.^2/b^2).>0.99
    ux[outsideellipse].=0.
    uy[outsideellipse].=0.
    uphi[outsideellipse].=0.

    extr = abs.(extrema(uphi[(X.^2/a^2+Y.^2/b^2).<.97]))
    cbi = findmax(extr)[2]
    maxval = (vmax==0) ? extr[cbi] : vmax

    streamplot(Float64.(X), Float64.(Y), Float64.(ux),
            Float64.(uy) ; color=Float64.(uphi),
            norm = PyPlot.matplotlib.colors.Normalize(vmin=-abs(maxval),vmax=abs(maxval)), kwargs...)

    ellipsex=a*cos.(range(0,stop=2π,length=100))
    ellipsey=b*sin.(range(0,stop=2π,length=100))
    plot(ellipsex,ellipsey,"k",linewidth=1.5)
    xlim([-1.1a,1.1a])

    axis("equal")
    axis("off");
end


function plot_velocity_meridional_x(b,c,v1; ngrid=50, kwargs...)

    Y, Z = grid(range(-b,stop=b,length=ngrid),range(-c,stop=c,length=ngrid))
    uy =  real.([v1[2](Mire.x=>0,Mire.y=>yt,Mire.z=>zt) for (yt,zt) in zip(Y,Z)])
    uz =  real.([v1[3](Mire.x=>0,Mire.y=>yt,Mire.z=>zt) for (yt,zt) in zip(Y,Z)])

    u = .√(uy.^2+uz.^2)

    outsideellipse=(Y.^2/b^2+Z.^2/c^2).>=.99
    uy[outsideellipse].=0.
    uz[outsideellipse].=0.
    u[outsideellipse].=0.

    extr = abs.(extrema(uz[(Z.^2/c^2+Y.^2/b^2).<.99]))
    cbi = findmax(extr)[2]

    streamplot(Float64.(Y), Float64.(Z), Float64.(uy),
            Float64.(uz) ; color=Float64.(uz), norm = PyPlot.matplotlib.colors.Normalize(vmin=-abs(extr[cbi]),vmax=abs(extr[cbi])), kwargs...)


    ellipsex=b*cos.(range(0,stop=2π,length=100))
    ellipsey=c*sin.(range(0,stop=2π,length=100))
    PyPlot.plot(ellipsex,ellipsey,"k",linewidth=2)
    PyPlot.axis("equal")
    PyPlot.axis("off");
end

function plot_velocity_meridional_y(a,c,v1; ngrid=50, kwargs...)

    X, Z = grid(range(-a,stop=a,length=ngrid),range(-c,stop=c,length=ngrid))
    ux =  real.([v1[1](Mire.x=>xt,Mire.y=>0,Mire.z=>zt) for (xt,zt) in zip(X,Z)])
    uz =  real.([v1[3](Mire.x=>xt,Mire.y=>0,Mire.z=>zt) for (xt,zt) in zip(X,Z)])

    u = .√(ux.^2+uz.^2)

    outsideellipse=(X.^2/a^2+Z.^2/c^2).>=.99
    ux[outsideellipse].=0.
    uz[outsideellipse].=0.
    u[outsideellipse].=0.

     extr = abs.(extrema(uz[(Z.^2/c^2+X.^2/a^2).<.99]))
    cbi = findmax(extr)[2]

    streamplot(Float64.(X), Float64.(Z), Float64.(ux),
            Float64.(uz) ; color=Float64.(uz), norm = PyPlot.matplotlib.colors.Normalize(vmin=-abs(extr[cbi]),vmax=abs(extr[cbi])), kwargs...)


    ellipsex=a*cos.(range(0,stop=2π,length=100))
    ellipsey=c*sin.(range(0,stop=2π,length=100))
    PyPlot.plot(ellipsex,ellipsey,"k",linewidth=2)
    PyPlot.axis("equal")
    PyPlot.axis("off");
end

df641 = one(Double64)

a,b,c,Le = df64"1.25",df64"0.8",df641,df64"1e-5"

b0f = (a,b,c)->b0_1_1(a,b,c)+b0_1_3(a,b,c)
pAform = (x^0*y^0+x)/df64"3"
pAform2 = y^0*x/df64"3"

b0Af = (a,b,c)-> b0_Aform(pAform,a,b,c);

mqg = ModelSetup(a,b,c,Le,b0Af, "ellipse_aform7", 7, QG())
mhyb = ModelSetup(a,b,c,Le,b0_2_8, "ellipse_b289", 9, Hybrid());


fname = joinpath(datapath,"QG_ellipse_aform7_df64_N7.jld")
JLD2.@load fname ω us

ω[1e-2.<imag.(ω).<10]

inds = reverse(eachindex(ω)[1e-2.<imag.(ω).<10])

PyPlot.rc("text",usetex=true)

inds=reverse(eachindex(ω)[1e-2.<imag.(ω).<10])
for i in 1:3
    figure(figsize=(3,2))
    plot_velocity_meridional_y(mqg.a,mqg.c,us[inds[i]]*1e5 .+eps(),cmap=:coolwarm,density=0.8);
    cb=colorbar();
    cb.set_label(L"u_z\, [10^{-5}]")
    if SAVEFIGS
        figname=joinpath(figpath,"streamplots_qg_meridional_$i.pdf")
        savefig(figname,bbox_inches="tight")
    end
end


inds=reverse(eachindex(ω)[1e-2.<imag.(ω).<10])
for i in 1:3
    f,ax=subplots(1,figsize=(3,2))
    plot_velocity_equator_uphi(mqg.a,mqg.b,us[inds[i]].+eps(),cbar=true,cmap=:coolwarm,density=0.8);
    cb=colorbar()
    cb.set_label(L"u_\varphi")
    if SAVEFIGS
        figname=joinpath(figpath,"streamplots_qg_equator_$i.pdf")
        savefig(figname,bbox_inches="tight")
    end

end

fname = joinpath(datapath,"Hybrid_ellipse_b289_df64_N9.jld")
JLD2.@load fname ω us;

ω[1e-2.<imag.(ω).<10]

inds=reverse(eachindex(ω)[1e-2.<imag.(ω).<10])
for i in 1:3
    figure(figsize=(3,2))
    plot_velocity_meridional_y(mhyb.a,mhyb.c,us[inds[i]]*1e5 .+eps(),cmap=:coolwarm,density=0.8);
    cb=colorbar();
    cb.set_label(L"u_z\, [10^{-5}]")
    if SAVEFIGS
        figname=joinpath(figpath,"streamplots_hyb_meridional_$i.pdf")
        savefig(figname,bbox_inches="tight")
    end
end


# f,ax=subplots(3,figsize=(4,8))
inds=reverse(eachindex(ω)[1e-2.<imag.(ω).<10])
for i in 1:3
    f,ax=subplots(1,figsize=(3,2))
    plot_velocity_equator_uphi(mhyb.a,mhyb.b,us[inds[i]].+eps(),cbar=true,cmap=:coolwarm,density=0.8);
    cb=colorbar()
    cb.set_label(L"u_\varphi")
    if SAVEFIGS
        figname=joinpath(figpath,"streamplots_hyb_equator_$i.pdf")
        savefig(figname,bbox_inches="tight")
    end

end
