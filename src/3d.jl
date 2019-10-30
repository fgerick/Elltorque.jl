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


function run(m::ModelSetup,datapath="",SAVEDATA=false)
    a,b,c,Le,b0 = m.a,m.b,m.c,m.Le,b.b0
    Ω = [0,0,1/Le]
    LHS,RHS,vs = Mire.assemblemhd(N,a,b,c,Ω,b0)
    A,B= complex.(Matrix(RHS)),complex.(Matrix(LHS))
    cmat = cacheint(N,a,b,c)
    S=GenericSchur.schur(A,B)
    ω=S.values
    evecs=eigvecs(S)
    us,bs=get_ub(evecs,vs,N,cmat)
    if SAVEDATA
        fname="/home/gerickf/data/eigenmodes_128bit_ell_N$(N).jld"
        JLD2.@save fname A B vs S ω evecs a b c Le Ω b0 us bs N
end
