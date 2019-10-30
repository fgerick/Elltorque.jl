using Revise, Test, Elltorque, Mire, DoubleFloats

#sphere float64
m1 = ModelSetup(1.,1.,1.,1e-5,b0_1_1)

#ellipsoid a,b,c=1.25,0.8,1.0
# big1=one(BigFloat)
# m2 = ModelSetup(big"1.25",big"0.8",big1,big"1e-5",(a,b,c)->b0_1_1(a,b,c)+b0_1_3(a,b,c))

df641 = one(Double64)
m2 = ModelSetup(df64"1.25",df64"0.8",df641,df64"1e-5",(a,b,c)->b0_1_1(a,b,c)+b0_1_3(a,b,c))



@testset "Run and save models" begin
    @test
    # Write your own tests here.
end
