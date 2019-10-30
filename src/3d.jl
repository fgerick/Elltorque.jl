function get_ub(evecs,vs,N,cmat)
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


function calculatemodes(m::ModelSetup,datapath="",SAVEDATA=false)
    a, b, c, Le, b0, N = m.a, m.b, m.c, m.Le, m.b0, m.N
    Ω = [0,0,1/Le]
    LHS, RHS, vs = Mire.assemblemhd(N, a, b, c, Ω, b0)
    A, B = complex.(Matrix(RHS)), complex.(Matrix(LHS))
    cmat = cacheint(N, a, b, c)
    S = GenericSchur.schur(A, B)
    ω = S.values
    evecs = eigvecs(S)
    us,bs = get_ub(evecs, vs, N, cmat)
    if SAVEDATA
        fname = joinpath(datapath,"3d_$(m.name)_$(typeof(a))_N$(N).jld")
        JLD2.@save fname A B vs S ω evecs m Ω us bs
    end
    return true
end
