
using Distributed
addprocs(parse(Int,ARGS[1]))
@everywhere using Elltorque, DoubleFloats
@show nprocs()
datapath = ARGS[2]
nle = parse(Int,ARGS[3])
lemax = parse(Double64,ARGS[4])
imagfield = parse(Int,ARGS[5])
modeldim = parse(Int,ARGS[6])
# Le = parse(Double64,ARGS[7])

@show nprocs()


@everywhere begin

    SAVEDATA = true
    df641 = one(Double64)

    #earth:
    # r0=df64"3480e3"
    # b = (r0+df64"1e4")/r0
    # a = 1/b
    # c = df641
    # Le = df64"0.0009320333592119371"

    a,b,c,Le = df64"1.25",df64"0.8",df641,df64"1e-5"

    b0f = (a,b,c)->b0_1_1(a,b,c)+b0_1_3(a,b,c)
    pAform = (x^0*y^0+x)/df64"3"
    b0Af = (a,b,c)-> b0_Aform(pAform,a,b,c)


    imdim = remotecall_fetch(()->modeldim,1)
    MDIM = (imdim == 1) ? QG() : ((imdim == 2) ? Hybrid() : Full())
    IMAG = remotecall_fetch(()->imagfield,1)

    if IMAG==1
        m0 = ModelSetup(df641,df641,df641,Le, b0_1_3,"malkussphere",3, Full())
    elseif IMAG==2
        m0 = ModelSetup(a,b,c,Le, b0_1_3,"malkusellipse",3, Full())
    elseif IMAG==3
        m0 = ModelSetup(a,b,c,Le,b0f, "ellipse1", 3, Full())
    elseif IMAG==4
        m0 = ModelSetup(a,b,c,Le,b0f, "ellipse2", 5, Full())
    elseif IMAG==5
        m0 = ModelSetup(a,b,c,Le,b0_2_6, "ellipse3", 5, Full())
    elseif IMAG==6
        m0 = ModelSetup(a,b,c,Le,b0Af, "ellipse4", 5, Full())
    elseif IMAG == 7
        b0f = b0Af
        m0 = ModelSetup(a,b,c,Le,b0f, "aform_ellipse2", 7, MDIM)
    elseif IMAG == 10
        b0f = b0_2_8
        m0 = ModelSetup(a,b,c,Le,b0f, "b028_ellipse", 9, MDIM)
    end

    T=Double64
    D=typeof(MDIM)
    dtypename="df64"
    Les=10.0.^(-10.0.^range(log10.(7),0.0,length=remotecall_fetch(()->nle, 1) ))
    # Les=10.0.^range(-7,-1,length=remotecall_fetch(()->nle,1))
    cmat = Mire.cacheint(m0.N,m0.a,m0.b,m0.c)
end

@time @sync @distributed for i=1:length(Les)

    m=ModelSetup{T,D}(m0.a,m0.b,m0.c,Les[i],m0.b0,"le_$i",m0.N)
    # Elltorque.calculatemodes(m,datapath,SAVEDATA,"df64")
    fname = joinpath(datapath,string(D)*"_$(m.name)_"*dtypename*"_N$(m.N).jld")
    if D == Full
        JLD2.@load fname A B vs S ω evecs m Ω us bs
    elseif D == Hybrid
        JLD2.@load fname A B vs vs_qg S ω evecs m Ω us bs
    elseif D == QG
        JLD2.@load fname A B vs_qg S ω evecs m Ω us bs
    end
    Ls = [Elltorque.angularmom(u,cmat,3) for u in us]
    fnameL = joinpath(datapath,"L_"*string(D)*"_$(m.name)_"*dtypename*"_N$(m.N).jld")
    JLD2.@save fnameL Ls
end
