function get_ub(evecs,vs,cmat)
    nev=size(evecs,2)
    us,bs=[],[]
    for i=1:nev
        # u = Mire.eigenvel(vs,evecs[1:end÷2,i]/mean(evecs[:,i]))
        # b = Mire.eigenvel(vs,evecs[end÷2+1:end,i]/mean(evecs[:,i]))
        u = Mire.eigenvel(vs,evecs[1:end÷2,i])#/mean(evecs[:,i]))
        b = Mire.eigenvel(vs,evecs[end÷2+1:end,i])#/mean(evecs[:,i]))

        energyu = Mire.inner_product(cmat,u,u)
        energyb = Mire.inner_product(cmat,b,b)
        energy = sqrt(energyu+energyb)
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
        u = Mire.eigenvel(vs_qg,evecs[1:nqg,i])#/mean(evecs[:,i]))
        b = Mire.eigenvel(vs,evecs[nqg+1:end,i])#/mean(evecs[:,i]))
        energyu = Mire.inner_product(cmat,u,u)
        energyb = Mire.inner_product(cmat,b,b)
        energy = sqrt(energyu+energyb)
        push!(us,u/energy)
        push!(bs,b/energy)
    end
    return us,bs
end


function tracking_b_3d(N,a,bs,c,α,Le,s0,s1,σ0,LHS0,RHS0,cmat,b0f; verbose=false, kwargs...)
    k=0
    bsout=[bs[1]]
    λs,us = [],[]
    target = σ0
    targets = [σ0]
    utarget = []
    LHS = copy(LHS0)
    RHS = copy(RHS0)
    vss = []

    db = bs[2]-bs[1]
    dbt = db
    bt = bs[1]
    iter = 0

    while bsout[end] + db >= bs[end]

        if (k == 0 )
            b = bs[1]
        else
            b = bsout[end]+db
        end

        if b < bt
            b = bt
        end

        a=1/b
        # cmat = Mire.cacheint(N,a,b,c)
        b0 = b0fun(a,b,c) #B0_malkus(a,b).+α*B0_x(a,b,c);

        LHS,RHS,vs = Mire.assemblemhd(N,a,b,c,1/Le * ez, b0)


        # λ,u = eigstarget(RHS, LHS, target; v0 = utarget, kwargs...)
        λ,u = eigstarget(Float64.(RHS), Float64.(LHS), target; v0=utarget, kwargs...)
        nev = length(λ)
        if k == 0
            imax = 1
            bt = bs[1] + db
        else
            corrs = abs.([cor(real.(u[:,i]/mean(u[:,i])),real.(utarget/mean(utarget))) for i=1:nev])
            corsort=sortperm(abs.(corrs),rev=true)
            cors=corrs[corsort[1:nev]]
            max_corr,imax = findmin(abs.(1.0 .- corrs))
            if corrs[imax] < 0.99 #if no correlating eigenvector is found the parameter stepping is too high
               if verbose
                    @warn "Correlation is only $(corrs[imax])!, lowering step"
                end
                dbt /= 2
                bt = bsout[end] + dbt
                flush(stdout)
                flush(stderr)
                continue
            else
                push!(bsout,b)
                bt = bsout[end] + db
                dbt = db
            end
        end

        push!(λs,λ[imax])
        push!(us,u[:,imax])
        push!(vss,vs)
        target = λ[imax]
        utarget= u[:,imax]
        if verbose
            @show b
            ω=abs(target)
            @show ω
            flush(stdout)
            flush(stderr)
        end
        push!(targets,target)
        k+=1
    end
    return λs,us,bsout,vss
end

function tracking_lehnert_3d(N,a,b,c,les,σ0,LHS0,RHS0,cmat, vs; verbose=false, kwargs...)
    k=0
    lesout=[les[1]]
    λs,us = [],[]
    target = σ0
    targets = [σ0]
    utarget = []
    LHS = deepcopy(LHS0)
    RHS = deepcopy(RHS0)

    dle = les[1] #diff(les) # les[2]-les[1]
    dlet = dle
    le_t = les[1]
    iter = 0
    mpeakold = 0
    while lesout[end] + dle <= les[end]

        if (k == 0 )
            Le = les[1]
        else
            Le = lesout[end]+dle #dle[k+1]
        end

        if Le>le_t
            Le = le_t
        end
        Ωvec = 1/Le*Mire.ez
        #only coriolis force depends on Le:
        RHS[1:end÷2,1:end÷2] = Mire.projectforce(N,cmat,vs,coriolis,Ωvec)
        λ,u = eigstarget(Float64.(RHS), Float64.(LHS), target; v0=utarget, kwargs...)
        nev = length(λ)
        if k == 0
            imax = 1
            le_t = les[1] + dle
        else
            corrs = abs.([cor(real.(u[:,i]/mean(u[:,i])),real.(utarget/mean(utarget))) for i=1:nev])
            corsort=sortperm(abs.(corrs),rev=true)
            cors=corrs[corsort]
            max_corr,imax = findmin(1.0 .- abs.(corrs))
            if abs(corrs[imax]) < 0.98 #if no correlating eigenvector is found the parameter stepping is too high
                if verbose
                    @warn "Correlation is only $(abs(corrs[imax]))!, lowering step"
                flush(stdout)
                flush(stderr)
                end
                dlet /= 2
                le_t = lesout[end] + dlet
                continue
            else
                push!(lesout,Le)
                dle = 2*(lesout[end]-lesout[end-1])
                le_t = lesout[end] + dle #dle[k+1]
                dlet = dle #[k+1]
            end
        end

        push!(λs,λ[imax])
        push!(us,u[:,imax])
        target = λ[imax]
        utarget= u[:,imax]

        if verbose
            ω=abs(target)
            @show ω

            flush(stdout)
            flush(stderr)
        end
        push!(targets,target)
        k+=1
        if verbose
            @show Le
            flush(stdout)
            flush(stderr)
        end
    end
    return λs,us,lesout
end

function eigstarget(A,B,target; kwargs...)
    P = lu(A-target*B)
    LO = LinearMap{ComplexF64}((y,x)->ldiv!(y,P,B*x),size(A,2))
    evals,u = eigs(LO; kwargs...)
    λ = 1 ./evals .+ target
    return λ,u
end

function eigen2(A,B; kwargs...)
    C = inv(Matrix(B))*A
    S = GenericSchur.triangularize(GenericSchur.gschur(C; kwargs...))
    u = eigvecs(S)
    return LinearAlgebra.Eigen(S.values, u)
end


function calculatemodes(m::ModelSetup{T,D},datapath="",SAVEDATA=false,dtypename="f64") where {T,D <: ModelDim}
    a, b, c, Le, b0, N = m.a, m.b, m.c, m.Le, m.b0, m.N
    Ω = [0,0,1/Le]

    if D==Full
        LHS, RHS, vs = Mire.assemblemhd(N, a, b, c, Ω, b0)
    elseif D==Hybrid
        LHS, RHS, vs, vs_qg = Mire.assemblemhd_hybrid(N, N, a, b, c, Ω, b0)
    elseif D==QG
        LHS, RHS, vs_qg = Mire.assemblemhd_qg(N, a, b, c, Ω, b0)
    else
        error("model must be one of: Full, Hybrid or QG!")
    end

    cmat = cacheint(N, a, b, c)

    if T <: LinearAlgebra.BlasFloat
        A,B = Matrix(RHS), Matrix(LHS)
        S = eigen(A, B)
    else
        A, B = complex.(Matrix(RHS)), complex.(Matrix(LHS))
        S = GenericSchur.schur(A, B)
    end

    ω = S.values
    evecs = eigvecs(S)

    if D == Full
        us, bs = get_ub(evecs, vs, cmat)
    elseif D == Hybrid
        us, bs = get_ub(evecs, vs, vs_qg, cmat)
    elseif D == QG
        us, bs = get_ub(evecs, vs_qg, cmat)
    end
    if SAVEDATA
        fname = joinpath(datapath,string(D)*"_$(m.name)_"*dtypename*"_N$(m.N).jld")
        if D == Full
            JLD2.@save fname A B vs S ω evecs m Ω us bs
        elseif D == Hybrid
            JLD2.@save fname A B vs vs_qg S ω evecs m Ω us bs
        elseif D == QG
            JLD2.@save fname A B vs_qg S ω evecs m Ω us bs
        end

    end
    return true
end

function calculatemodesthread(mutex,m::ModelSetup{T,D},datapath="",SAVEDATA=false,dtypename="f64", torques=true) where {T,D <: ModelDim}
    a, b, c, Le, b0, N = m.a, m.b, m.c, m.Le, m.b0, m.N
    Ω = [0,0,1/Le]

    if D==Full
        LHS, RHS, vs = Mire.assemblemhd(N, a, b, c, Ω, b0)
    elseif D==Hybrid
        LHS, RHS, vs, vs_qg = Mire.assemblemhd_hybrid(N, N, a, b, c, Ω, b0)
    elseif D==QG
        LHS, RHS, vs_qg = Mire.assemblemhd_qg(N, a, b, c, Ω, b0)
    else
        error("model must be one of: Full, Hybrid or QG!")
    end

    cmat = cacheint(N, a, b, c)

    if T <: LinearAlgebra.BlasFloat
        A,B = Matrix(RHS), Matrix(LHS)
        S = eigen(A, B)
    else
        A, B = complex.(Matrix(RHS)), complex.(Matrix(LHS))
        S = GenericSchur.schur(A, B)
    end

    ω = S.values
    evecs = eigvecs(S)

    if D == Full
        us, bs = get_ub(evecs, vs, cmat)
    elseif D == Hybrid
        us, bs = get_ub(evecs, vs, vs_qg, cmat)
    elseif D == QG
        us, bs = get_ub(evecs, vs_qg, cmat)
    end


    ugua = split_ug_ua.(us,m.a,m.b,m.c)
    ug,ua = getindex.(ugua,1), getindex.(ugua,2)
    # end
    # Γp, Γptot, Γpmag, Lω, Lωa, Γem, Γcor, Γpageo = torquebalance(m.N,m.a,m.b,m.c,us,ug,ua,bs,ω,m.b0,Ω)
    if torques
        Γp, Γptot, Γpmag, Lω, Γem, Γcor = torquebalance(m.N,m.a,m.b,m.c,us,ug,ua,bs,ω,m.b0,Ω)
    end
    # Γp_uni, Γptot_uni, Γpmag_uni, Lω_uni, Γem_uni, Γcor_uni = torquebalance(m.N,m.a,m.b,m.c,us,ug,ua,bs,ω,m.b0,Ω,takeuni=true)



    lock(mutex)
    if SAVEDATA
        fname = joinpath(datapath,string(D)*"_$(m.name)_"*dtypename*"_N$(m.N).jld")
        if D == Full
            JLD2.@save fname A B vs S ω evecs m Ω us bs
        elseif D == Hybrid
            JLD2.@save fname A B vs vs_qg S ω evecs m Ω us bs
        elseif D == QG
            JLD2.@save fname A B vs_qg S ω evecs m Ω us bs
        end
    end

    if SAVEDATA && torques
        fname = joinpath(datapath,"torquebalance_"*string(D)*"_$(m.name)_"*dtypename*"_N$(m.N).jld")
        JLD2.@save fname Γp Γptot Γpmag Lω Γem Γcor
        # fname = joinpath(datapath,"torquebalance_uni_"*string(D)*"_$(m.name)_"*dtypename*"_N$(m.N).jld")
        # JLD2.@save fname Γp_uni Γptot_uni Γpmag_uni Lω_uni Γem_uni Γcor_uni
    end
    unlock(mutex)
    return true
end



function runcalculations(SAVEDATA,datapath)
    ## QG models

    df641 = one(Double64)

    #earth:
    # r0=df64"3480e3"
    # b = (r0+df64"1e4")/r0
    # a = 1/b
    # c = df641
    # Le = df64"0.0009320333592119371"
    a,b,c,Le = df64"1.25",df64"0.8",df641,df64"1e-5"
    # a,b,c,Le = df641,df641,df641,df64"1e-5"

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

    # for m in mall
    #     calculatemodes(m,datapath,SAVEDATA,"df64")
    #     loadandcalculatetorque(m,datapath,SAVEDATA,"df64")
    # end

    mutex = Threads.Mutex()
    Threads.@threads for m in mall
        calculatemodesthread(mutex,m,datapath,SAVEDATA,"df64")
    end


end





function runcalculationslehnert(m0::ModelSetup{T,D},Les,SAVEDATA,datapath) where {T,D<:ModelDim}


    mutex = Threads.Mutex()
    Threads.@threads for i=1:length(Les)
        m=ModelSetup{T,D}(m0.a,m0.b,m0.c,Les[i],m0.b0,"le_$i",m0.N)
        calculatemodesthread(mutex,m,datapath,SAVEDATA,"df64",false)
    end


end
