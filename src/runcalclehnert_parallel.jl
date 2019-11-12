
using Distributed
addprocs(parse(Int,ARGS[1]))
datapath = ARGS[2]
nle = parse(Int,ARGS[3])
imagfield = parse(Int,ARGS[4])

@show nprocs()
@everywhere using Elltorque, DoubleFloats


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


    ## 3D models

    ## 3D models
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
    end

    T=Double64
    D=Full
    Les=10.0.^(-10.0.^range(log10.(7),0.0,length=remotecall_fetch(()->nle))
    # Les=10.0.^range(-7,-1,length=remotecall_fetch(()->nle,1))
end

@time @sync @distributed for i=1:length(Les)

    m=ModelSetup{T,D}(m0.a,m0.b,m0.c,Les[i],m0.b0,"le_$i",m0.N)
    Elltorque.calculatemodes(m,datapath,SAVEDATA,"df64")
end
