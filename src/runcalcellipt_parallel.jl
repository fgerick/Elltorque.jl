
using Distributed
addprocs(parse(Int,ARGS[1]))
@everywhere using Elltorque, DoubleFloats
@show nprocs()
datapath = ARGS[2]
nb = parse(Int,ARGS[3])
bmax = parse(Double64,ARGS[4])
imagfield = parse(Int,ARGS[5])
modeldim = parse(Int,ARGS[6])

@everywhere begin

    SAVEDATA = true
    df641 = one(Double64)

    #earth:
    # r0=df64"3480e3"
    # b = (r0+df64"1e4")/r0
    # a = 1/b
    # c = df641
    # Le = df64"0.0009320333592119371"

    a,b,c,Le = df641,df641,df641,df64"1e-4"

    b0f = (a,b,c)->b0_1_1(a,b,c)+b0_1_3(a,b,c)
    pAform = (x^0*y^0+x)/df64"3"
    b0Af = (a,b,c)-> b0_Aform(pAform,a,b,c)

    ## 3D models
    imdim = remotecall_fetch(()->madeldim,1)
    MDIM = (imdim == 1) ? QG() : ((imdim == 2) ? Hybrid() : Full())
    IMAG = remotecall_fetch(()->imagfield,1)
    if IMAG == 1
        m0 = ModelSetup(df641,df641,df641,Le, b0_1_3,"malkussphere", 3, MDIM)
    elseif IMAG == 2
        m0 = ModelSetup(a,b,c,Le, b0_1_3, "malkusellipse", 3, MDIM)
    elseif IMAG == 3
        m0 = ModelSetup(a,b,c,Le,b0f, "ellipse1", 3, MDIM)
    elseif IMAG == 4
        m0 = ModelSetup(a,b,c,Le,b0f, "ellipse2", 5, MDIM)
    elseif IMAG == 5
        m0 = ModelSetup(a,b,c,Le,b0_2_6, "b026_ellipse", 5, MDIM)
    elseif IMAG == 6
        m0 = ModelSetup(a,b,c,Le,b0Af, "aform_ellipse", 5, MDIM)
    elseif IMAG == 7
        m0 = ModelSetup(a,b,c,Le,b0Af, "aform_ellipse2", 7, MDIM)
    elseif IMAG == 8
        m0 = ModelSetup(a,b,c,Le,b0_2_6, "b026_ellipse2", 7, MDIM)
    end

    T = Double64
    D = Full
    # bs = df64"10.0".^range(0,log10.(remotecall_fetch(()->bmax,1)),length=remotecall_fetch(()->nb,1))
    # bs = range(df641,remotecall_fetch(()->bmax,1),length=remotecall_fetch(()->nb,1))
    ϵs = df64"10.0".^range(-7,log10.(remotecall_fetch(()->bmax,1)),length=remotecall_fetch(()->nb,1))
end

@time @sync @distributed for i=1:length(ϵs)
    ϵ = ϵs[i]
    b = ((1 - ϵ)/(1 + ϵ))^(1//4)
    a = 1/b
    b0 = b0Af(a,b,m0.c)
    m=ModelSetup{T,D}(a,b,m0.c,m0.Le,b0,"eps_$i",m0.N)
    Elltorque.calculatemodes(m,datapath,SAVEDATA,"df64")
end
