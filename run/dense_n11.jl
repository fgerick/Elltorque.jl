#calculates the data for Figure 12
using Distributed
addprocs(parse(Int,ARGS[1]))
@show nprocs()
@everywhere using DoubleFloats, Elltorque, JLD2
datapath = ARGS[2]


@everywhere begin

    SAVEDATA = true
    nle = 200
    df641 = one(Double64)

    a,b,c,Le = df64"1.25",df64"0.8",df641,df64"1e-7"

    pAform = (x^0*y^0+x)/df64"3"
    b0Af = (a,b,c)-> b0_Aform(pAform,a,b,c)

    m0 = ModelSetup(a,b,c,Le,b0Af, "aform_ellipse_n11",11,QG() )

    T = Double64
    D = QG
    dtypename = "df64"
    Les = df64"10.0".^range(-7,0,length=nle)
    cmat = Mire.cacheint(m0.N,m0.a,m0.b,m0.c)
end

@time @sync @distributed for i=1:length(Les)

    m = ModelSetup{T,D}(m0.a,m0.b,m0.c,Les[i],m0.b0,"le_$i",m0.N)
    Elltorque.calculatemodes(m,datapath,SAVEDATA,dtypename)
    fname = joinpath(datapath,string(D)*"_$(m.name)_"*dtypename*"_N$(m.N).jld")
    JLD2.@load fname A B vs_qg S ω evecs m Ω us bs
    Ls = [Elltorque.angularmom(u,cmat,3) for u in us]
    fnameL = joinpath(datapath,"L_"*string(D)*"_$(m.name)_"*dtypename*"_N$(m.N).jld")
    JLD2.@save fnameL Ls Les
end
