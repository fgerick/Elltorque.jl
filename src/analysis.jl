function get_ub1(λs,vs,vs_qg,cmat; ekin=false, getenergies=false)
    nqg=length(vs_qg)
    u = Mire.eigenvel(vs_qg,λs[1:nqg])#/mean(evecs[:,i]))
    b = Mire.eigenvel(vs,λs[nqg+1:end])#/mean(evecs[:,i]))
    energyu = Mire.inner_product(cmat,u,u)
    energyb = Mire.inner_product(cmat,b,b)
    if ekin
        u = u/√energyu
        b = b/√energyb
    else
        en = sqrt(energyu + energyb)
        u = u/en
        b = b/en
    end
    if getenergies
        return u,b,energyu,energyb
    else
        return u,b
    end
end

function get_ub(evecs,vs,vs_qg,cmat; ekin=false, getenergies=false)
    nev = size(evecs,2)
    nqg = length(vs_qg)
    us,bs=[],[]
    if getenergies
        eks,ebs = [],[]
    end
    for i=1:nev
        if getenergies
            u,b,ek,eb = get_ub1(evecs[:,i],vs,vs_qg,cmat; ekin=ekin, getenergies=true)
            push!(eks,ek)
            push!(ebs,eb)
        else
            u,b = get_ub1(evecs[:,i],vs,vs_qg,cmat; ekin=ekin, getenergies=false)
        end
        push!(us,u)
        push!(bs,b)

    end

    if getenergies
        return us,bs,eks,ebs
    else
        return us,bs
    end
end

get_ub(evecs,vs,cmat; kwargs...) = get_ub(evecs,vs,vs,cmat; kwargs...)

function tracking_ellipt(m::ModelSetup{T,D},ϵs,σ0,LHS0,RHS0,b0f; verbose=false, kwargs...) where {T<:Number,D<:ModelDim}
    N,a,c,Le = m.N,m.a,m.c,m.Le
    Ω = [0,0,1/Le]
    k = 0
    ϵsout = [ϵs[1]]
    λs,us = [],[]
    target = σ0
    targets = [σ0]
    utarget = [] #initial guess eigenvector
    LHS = copy(LHS0)
    RHS = copy(RHS0)
    vss = []

    dϵ = ϵs[2]-ϵs[1]
    dϵt = dϵ
    ϵt = ϵs[1]
    iter = 0

    while ϵsout[end] < ϵs[end]-eps()

        if (k == 0 )
            ϵ = ϵs[1]
        else
            ϵ = ϵsout[end]+dϵ
        end

        if ϵ > ϵt
            ϵ = ϵt
        end

        if ϵ >= ϵs[end]
            ϵ = ϵs[end]
        end
        b = ((1 - ϵ)/(1 + ϵ))^(1//4)
        a = 1/b
        b0 = b0f(a,b,c) #B0_malkus(a,b).+α*B0_x(a,b,c);
        if D==Full
            LHS, RHS, vs = Mire.assemblemhd(N, a, b, c, Ω, b0)
        elseif D==Hybrid
            LHS, RHS, vs, vs_qg = Mire.assemblemhd_hybrid(N, N, a, b, c, Ω, b0)
        elseif D==QG
            LHS, RHS, vs = Mire.assemblemhd_qg(N, a, b, c, Ω, b0)
        else
            error("model must be one of: Full, Hybrid or QG!")
        end



        λ,u = eigstarget(Float64.(RHS), Float64.(LHS), target; v0=utarget, kwargs...)
        nev = length(λ)

        # find correlating eigenvector:
        if k == 0
            imax = 1
            ϵt = ϵs[1] + dϵ
        else
            corrs = abs.([cor(real.(u[:,i]/mean(u[:,i])),real.(utarget/mean(utarget))) for i=1:nev])
            corsort=sortperm(abs.(corrs),rev=true)
            cors=corrs[corsort[1:nev]]
            max_corr,imax = findmin(abs.(1.0 .- corrs))
            if corrs[imax] < 0.99 #if no correlating eigenvector is found the parameter stepping is too high
               if verbose
                    @warn "Correlation is only $(corrs[imax])!, lowering step"
                end
                dϵt /= 2
                ϵt = ϵsout[end] + dϵt
                flush(stdout)
                flush(stderr)
                continue
            else
                push!(ϵsout,ϵ)
                dϵ = 1.2*(ϵsout[end]-ϵsout[end-1])
                ϵt = ϵsout[end] + dϵ
                dϵt = dϵ
            end
        end

        push!(λs,λ[imax])
        push!(us,u[:,imax])
        if D==Full || D==QG
            push!(vss,vs)
        else
            push!(vss,(vs,vs_qg))
        end
        target = λ[imax]
        utarget = u[:,imax]
        if verbose
            @show ϵ
            ω = abs(target)
            @show ω
            # flush(stdout)
            # flush(stderr)
        end
        push!(targets,target)
        k+=1
    end
    return λs,us,ϵsout,vss
end

function tracking_ellipt_reverse(m::ModelSetup{T,D},ϵ0,dϵ,σ0,LHS0,RHS0,b0f; zerothresh=eps(), verbose=false, kwargs...) where {T<:Number,D<:ModelDim}
    N,a,c,Le = m.N,m.a,m.c,m.Le
    Ω = [0,0,1/Le]
    k = 0
    ϵsout = [ϵ0]
    λs,us = [],[]
    target = σ0
    targets = [σ0]
    utarget = [] #initial guess eigenvector
    LHS = copy(LHS0)
    RHS = copy(RHS0)
    vss,cmats = [],[]

    dϵt = dϵ
    ϵt = ϵ0
    iter = 0

    while ϵsout[end] > zerothresh

        if (k == 0 )
            ϵ = ϵ0
        else
            ϵ = ϵsout[end]-dϵ
        end

        if ϵ < ϵt
            ϵ = ϵt
        end

        if ϵ <= 0.0
            ϵ = 0.0
        end

        b = ((1 - ϵ)/(1 + ϵ))^(1//4)
        a = 1/b
        b0 = b0f(a,b,c) #B0_malkus(a,b).+α*B0_x(a,b,c);
        cmat = cacheint(N,a,b,c)
        if D==Full
            LHS, RHS, vs = Mire.assemblemhd(N, a, b, c, Ω, b0, cmat=cmat)
        elseif D==Hybrid
            LHS, RHS, vs, vs_qg = Mire.assemblemhd_hybrid(N, N, a, b, c, Ω, b0, cmat=cmat)
        elseif D==QG
            LHS, RHS, vs = Mire.assemblemhd_qg(N, a, b, c, Ω, b0, cmat=cmat)
        else
            error("model must be one of: Full, Hybrid or QG!")
        end



        λ,u = eigstarget(Float64.(RHS), Float64.(LHS), target; v0=utarget, kwargs...)
        nev = length(λ)

        # find correlating eigenvector:
        if k == 0
            imax = 1
            ϵt = ϵ0 - dϵ
        else
            corrs = abs.([cor(real.(u[:,i]/mean(u[:,i])),real.(utarget/mean(utarget))) for i=1:nev])
            corsort=sortperm(abs.(corrs),rev=true)
            cors=corrs[corsort[1:nev]]
            max_corr,imax = findmin(abs.(1.0 .- corrs))
            if corrs[imax] < 0.99 #if no correlating eigenvector is found the parameter stepping is too high
               if verbose
                    @warn "Correlation is only $(corrs[imax])!, lowering step"
                end
                dϵt /= 2
                ϵt = ϵsout[end] - dϵt
                flush(stdout)
                flush(stderr)
                continue
            else
                push!(ϵsout,ϵ)
                dϵ = 1.2*abs(ϵsout[end-1]-ϵsout[end])
                dϵ = (dϵ > ϵsout[end]/2) ? ϵsout[end]/2 : dϵ
                ϵt = ϵsout[end] - dϵ
                dϵt = dϵ
            end
        end

        push!(λs,λ[imax])
        push!(us,u[:,imax])
        push!(cmats,cmat)
        if D==Full || D==QG
            push!(vss,vs)
        else
            push!(vss,(vs,vs_qg))
        end
        target = λ[imax]
        utarget = u[:,imax]
        if verbose
            @show ϵ
            ω = abs(target)
            @show ω
            # flush(stdout)
            # flush(stderr)
        end
        push!(targets,target)
        k+=1
    end
    return λs,us,ϵsout,vss,cmats
end

function tracking_lehnert(m::ModelSetup{T,D},les,σ0,LHS0,RHS0; verbose=false, corrtol=0.99, maxdlefac=100, kwargs...) where {T<:Number,D<:ModelDim}
    k=0
    lesout=[les[1]]
    λs,us = [],[]
    target = σ0
    targets = [σ0]
    utarget = []
    LHS = copy(LHS0)
    RHS = copy(RHS0)

    dle = les[1]/maxdlefac #diff(les) # les[2]-les[1]
    dlet = dle
    le_t = les[1]
    iter = 0
    mpeakold = 0
    a,b,c = m.a,m.b,m.c
    N = m.N
    cmat = cacheint(N,a,b,c)
    if D==Full
        vs = vel(N,a,b,c)
    else
        vs = Mire.qg_vel(N,a,b,c)
    end

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
        Mire.projectforce!(view(RHS,1:length(vs),1:length(vs)),cmat,vs,vs,coriolis,Ωvec)

        λ,u = Elltorque.eigstarget(Float64.(RHS), Float64.(LHS), target; v0=utarget, kwargs...)
        nev = length(λ)
        if k == 0
            imax = 1
            le_t = les[1] + dle
        else
            corrs = [abs(cor(u[:,i],utarget)) for i=1:nev]
            max_corr,imax = findmax(corrs)
            if abs(corrs[imax]) < corrtol #if no correlating eigenvector is found the parameter stepping is too high
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
                dle = 2*abs(lesout[end]-lesout[end-1])
                if dle > lesout[end]/maxdlefac
                    dle = lesout[end]/maxdlefac
                end
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

            # flush(stdout)
            # flush(stderr)
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


function calculatemodes(m::ModelSetup{T,D},datapath="",SAVEDATA=false,dtypename="f64";ekin=false) where {T,D <: ModelDim}
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
        us, bs = get_ub(evecs, vs, cmat, ekin=ekin)
    elseif D == Hybrid
        us, bs = get_ub(evecs, vs, vs_qg, cmat, ekin=ekin)
    elseif D == QG
        us, bs = get_ub(evecs, vs_qg, cmat, ekin=ekin)
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
