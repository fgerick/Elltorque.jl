using Revise, Test, Elltorque, Mire

m1 = ModelSetup(1.,1.,1.,1e-5,b0_1_1)

big1=one(BigFloat)
m2 = ModelSetup(big1,big1,big1,big"1e-5",(a,b,c)->b0_1_1(a,b,c)+b0_1_3(a,b,c))
