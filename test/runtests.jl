using Elltorque, Mire, DoubleFloats
using Test



SAVEDATA=true
datapath=""






## 3D models

df641 = one(Double64)
a,b,c,Le = df64"1.25",df64"0.8",df641,df64"1e-5"
b0f = (a,b,c)->b0_1_1(a,b,c)+b0_1_3(a,b,c)
m1 = ModelSetup(df641,df641,df641,Le, b0_1_3,"malkussphere",3, Full())
m2 = ModelSetup(a,b,c,Le, b0_1_3,"malkusellipse",3, Full())
m3 = ModelSetup(a,b,c,Le,b0f, "ellipse1", 3, Full())
# m4 = ModelSetup(a,b,c,Le,b0f, "ellipse2", 5, Full())
# m5 = ModelSetup(a,b,c,Le,b0_2_6, "ellipse3", 5, Full())
# m6 = ModelSetup(a,b,c,Le,b0_Aform, "ellipse3", 5, Full())

@testset "3D models" begin
    for m in [m1,m2,m3] #,m4,m5]
        @test calculatemodes(m,datapath,SAVEDATA)
        @test loadandcalculatetorque(m,datapath,SAVEDATA)
    end
end



## Hybrid models

df641 = one(Double64)
a,b,c,Le = df64"1.25",df64"0.8",df641,df64"1e-5"
b0f = (a,b,c)->b0_1_1(a,b,c)+b0_1_3(a,b,c)
m1 = ModelSetup(df641,df641,df641,Le, b0_1_3,"malkussphere",3, Hybrid())
m2 = ModelSetup(a,b,c,Le, b0_1_3,"malkusellipse",3, Hybrid())
m3 = ModelSetup(a,b,c,Le,b0f, "ellipse1", 3, Hybrid())
# m4 = ModelSetup(a,b,c,Le,b0f, "ellipse2", 5, Hybrid())
# m5 = ModelSetup(a,b,c,Le,b0_2_6, "ellipse3", 5, Hybrid())
# m6 = ModelSetup(a,b,c,Le,b0_Aform, "ellipse3", 5, Full())

@testset "Hybrid models" begin
    for m in [m1,m2,m3] #,m4,m5]
        @test calculatemodes(m,datapath,SAVEDATA)
        @test loadandcalculatetorque(m,datapath,SAVEDATA)
    end
end
