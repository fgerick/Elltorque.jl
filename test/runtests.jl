using Elltorque, Mire, DoubleFloats
using Test


m1 = ModelSetup(1.,1.,1.,1e-5,b0_1_1,"sphere",3, Full())

#ellipsoid a,b,c=1.25,0.8,1.0

df641 = one(Double64)
m2 = ModelSetup(df64"1.25",df64"0.8",df641,df64"1e-5",
                (a,b,c)->b0_1_1(a,b,c)+b0_1_3(a,b,c),
                "model2", 3, Full())

m3 = ModelSetup(df64"1.25",df64"0.8",df641,df64"1e-5",
                (a,b,c)->b0_1_1(a,b,c)+b0_1_3(a,b,c),
                "model2", 3, Hybrid())


SAVEDATA=true
datapath=""
@testset "Run and save models" begin
    @test calculatemodes(m1,datapath,SAVEDATA)
    @test calculatemodes(m2,datapath,SAVEDATA)
    @test calculatemodes(m3,datapath,SAVEDATA)

    # Write your own tests here.
end
