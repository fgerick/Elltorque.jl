function get_ub(evecs,vs,cmat)
    nev=size(evecs,2)
    us,bs=[],[]
    for i=1:nev
        u = Mire.eigenvel(vs,evecs[1:end÷2,i]/mean(evecs[:,i]))
        b = Mire.eigenvel(vs,evecs[end÷2+1:end,i]/mean(evecs[:,i]))
        energyu = √Mire.inner_product(cmat,u,u)
        energyb = √Mire.inner_product(cmat,b,b)
        energy = energyu+energyb
        push!(us,u/energy)
        push!(bs,b/energy)
    end
    return us,bs
end

function get_ub(evecs,vs,vs_qg,cmat)
    nev=size(evecs,2)
    nqg=length(vs_qg)
    us,bs=[],[]
    for i=1:nev
        u = Mire.eigenvel(vs_qg,evecs[1:nqg,i]/mean(evecs[:,i]))
        b = Mire.eigenvel(vs,evecs[nqg+1:end,i]/mean(evecs[:,i]))
        energyu = √Mire.inner_product(cmat,u,u)
        energyb = √Mire.inner_product(cmat,b,b)
        energy = energyu+energyb
        push!(us,u/energy)
        push!(bs,b/energy)
    end
    return us,bs
end


function calculatemodes(m::ModelSetup{T,D},datapath="",SAVEDATA=false) where {T,D <: ModelDim}
    a, b, c, Le, b0, N = m.a, m.b, m.c, m.Le, m.b0, m.N
    Ω = [0,0,1/Le]

    if D==Full
        LHS, RHS, vs = Mire.assemblemhd(N, a, b, c, Ω, b0)
    elseif D==Hybrid
        LHS, RHS, vs, vs_qg = Mire.assemblemhd_hybrid(N, N, a, b, c, Ω, b0)
    end

    A, B = complex.(Matrix(RHS)), complex.(Matrix(LHS))
    cmat = cacheint(N, a, b, c)
    S = GenericSchur.schur(A, B)
    ω = S.values
    evecs = eigvecs(S)

    if D == Full
        us, bs = get_ub(evecs, vs, cmat)
    elseif D == Hybrid
        us, bs = get_ub(evecs,vs,vs_qg, cmat)
    end
    if SAVEDATA
        fname = joinpath(datapath,string(D)*"_$(m.name)_$(T)_N$(m.N).jld")
        JLD2.@save fname A B vs S ω evecs m Ω us bs
    end
    return true
end



function runcalculations(SAVEDATA,datapath)

    ## 3D models

    df641 = one(Double64)
    a,b,c,Le = df64"1.25",df64"0.8",df641,df64"1e-5"
    b0f = (a,b,c)->b0_1_1(a,b,c)+b0_1_3(a,b,c)
    m1 = ModelSetup(df641,df641,df641,Le, b0_1_3,"malkussphere",3, Full())
    m2 = ModelSetup(a,b,c,Le, b0_1_3,"malkusellipse",3, Full())
    m3 = ModelSetup(a,b,c,Le,b0f, "ellipse1", 3, Full())
    m4 = ModelSetup(a,b,c,Le,b0f, "ellipse2", 5, Full())
    m5 = ModelSetup(a,b,c,Le,b0_2_6, "ellipse3", 5, Full())
    # m6 = ModelSetup(a,b,c,Le,b0_Aform, "ellipse3", 5, Full())

    for m in [m1,m2,m3,m4,m5]
         calculatemodes(m,datapath,SAVEDATA)
         loadandcalculatetorque(m,datapath,SAVEDATA)
    end



    ## Hybrid models

    df641 = one(Double64)
    a,b,c,Le = df64"1.25",df64"0.8",df641,df64"1e-5"
    b0f = (a,b,c)->b0_1_1(a,b,c)+b0_1_3(a,b,c)
    m1 = ModelSetup(df641,df641,df641,Le, b0_1_3,"malkussphere",3, Hybrid())
    m2 = ModelSetup(a,b,c,Le, b0_1_3,"malkusellipse",3, Hybrid())
    m3 = ModelSetup(a,b,c,Le,b0f, "ellipse1", 3, Hybrid())
    m4 = ModelSetup(a,b,c,Le,b0f, "ellipse2", 5, Hybrid())
    m5 = ModelSetup(a,b,c,Le,b0_2_6, "ellipse3", 5, Hybrid())
    # m6 = ModelSetup(a,b,c,Le,b0_Aform, "ellipse3", 5, Full())

    for m in [m1,m2,m3,m4,m5]
        calculatemodes(m,datapath,SAVEDATA)
        loadandcalculatetorque(m,datapath,SAVEDATA)
    end
end
