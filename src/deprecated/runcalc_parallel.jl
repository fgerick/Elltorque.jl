
using Distributed
addprocs(parse(Int,ARGS[1]))
datapath = ARGS[2]
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

    m1qg = ModelSetup(df641,df641,df641,Le, b0_1_3,"malkussphere",3, QG())
    m2qg = ModelSetup(a,b,c,Le, b0_1_3,"malkusellipse",3, QG())
    m3qg = ModelSetup(a,b,c,Le,b0Af, "ellipse_aform", 5, QG())

    msqg = [m1qg,m2qg,m3qg]

    ## Hybrid models

    m1h = ModelSetup(df641,df641,df641,Le, b0_1_3,"malkussphere",3, Hybrid())
    m2h = ModelSetup(a,b,c,Le, b0_1_3,"malkusellipse",3, Hybrid())
    m3h = ModelSetup(a,b,c,Le,b0f, "ellipse1", 3, Hybrid())
    m4h = ModelSetup(a,b,c,Le,b0f, "ellipse2", 5, Hybrid())
    m5h = ModelSetup(a,b,c,Le,b0_2_6, "ellipse3", 5, Hybrid())
    m6h = ModelSetup(a,b,c,Le,b0Af, "ellipse4", 5, Hybrid())

    mshybrid = [m1h,m2h,m3h,m4h,m5h,m6h]

    ## 3D models

    m1 = ModelSetup(df641,df641,df641,Le, b0_1_3,"malkussphere",3, Full())
    m2 = ModelSetup(a,b,c,Le, b0_1_3,"malkusellipse",3, Full())
    m3 = ModelSetup(a,b,c,Le,b0f, "ellipse1", 3, Full())
    m4 = ModelSetup(a,b,c,Le,b0f, "ellipse2", 5, Full())
    m5 = ModelSetup(a,b,c,Le,b0_2_6, "ellipse3", 5, Full())
    m6 = ModelSetup(a,b,c,Le,b0Af, "ellipse4", 5, Full())

    msfull = [m1,m2,m3,m4,m5,m6]

    mall = vcat(msqg,mshybrid,msfull)
end

@time @sync @distributed for m in mall
    calculatemodes(m,datapath,SAVEDATA,"df64")
    loadandcalculatetorque(m,datapath,SAVEDATA,"df64")
end
