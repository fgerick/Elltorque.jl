using Distributed
addprocs(parse(Int,ARGS[1]))
@everywhere using Elltorque, DoubleFloats, JLD2
@show nprocs()
datapath = ARGS[2]
nle = parse(Int,ARGS[3])
imagfield = parse(Int,ARGS[4])
modeldim = parse(Int,ARGS[5])
truncdegree = parse(Int,ARGS[6])


#earth:
# r0=df64"3480e3"
# b = (r0+df64"1e4")/r0
# a = 1/b
# c = df641
# Le = df64"0.0009320333592119371"
@everywhere begin

    SAVEDATA = true
    df641 = one(Double64)
    a,b,c,Le = df64"1.25",df64"0.8",df641,df64"1e-5"

    b0f = (a,b,c)->b0_1_1(a,b,c)+b0_1_3(a,b,c)
    pAform = (x^0*y^0+x)/df64"3"
    b0Af = (a,b,c)-> b0_Aform(pAform,a,b,c)


    imdim = remotecall_fetch(()->modeldim,1)
    MDIM = (imdim == 1) ? QG() : ((imdim == 2) ? Hybrid() : Full())
    IMAG = remotecall_fetch(()->imagfield,1)
    N = remotecall_fetch(()->truncdegree,1)

    if IMAG==1
        b0f = b0Af
        m0 = ModelSetup(a,b,c,Le,b0f, "aform_ellipse", N, MDIM)
    elseif IMAG == 2
        b0f = b0_2_8
        m0 = ModelSetup(a,b,c,Le,b0f, "b028_ellipse", N, MDIM)
    end

    T=Double64
    D=typeof(MDIM)
    dtypename="f64"
    les=10.0.^range(-4,log10.(5e-2),length=remotecall_fetch(()->nle, 1) )
    cmat = Mire.cacheint(m0.N,m0.a,m0.b,m0.c)
    m = m0
    N,a,b,c,b0 = m.N,m.a,m.b,m.c,m.b0
    if D == Full
        LHS, RHS, vs1 = Mire.assemblemhd(N, a, b, c, [0,0,1/m.Le], b0)
        vs2 = vs1
        elseif D == Hybrid
            LHS, RHS, vs1, vs_qg = Mire.assemblemhd_hybrid(N, N, a, b, c, [0,0,1/m.Le], b0)
        elseif D == QG
            LHS, RHS, vs1 = Mire.assemblemhd_qg(N, a, b, c, [0,0,1/m.Le], b0)
            vs2 = vs1
        else
            error("model must be one of: Full, Hybrid or QG!")
    end

    RHSt = copy(RHS)
    cmat = cacheint(N,a,b,c)
    invL = inv(Matrix(LHS))
    ωs, us, bs, eks, ebs = [], [], [], [], []
    @sync @distributed for i in 1:length(les)
        Mire.projectforce!(view(RHSt,1:length(vs1),1:length(vs2)),cmat,vs1,vs2,coriolis,[0,0,1/les[i]])
        esol = eigen(Float64.(invL*RHSt))

        us,bs,eks,ebs = Elltorque.get_ub(esol.vectors,vs1,vs2,Float64.(cmat); ekin=false, getenergies=true)
        Ls = [Elltorque.angularmom(u,cmat,3) for u in us]
        fnameL = joinpath(datapath,"Lle_$(i)_"*string(D)*"_$(m.name)_"*dtypename*"_N$(m.N).jld")
        JLD2.@save fnameL Ls les us bs eks ebs
    end

end
