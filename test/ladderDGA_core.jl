@testset "index(1:) to freq" begin
    #no shift tests
    @test all(BSE_SC.OneToIndex_to_Freq(1,1,5,5,0) .== (-5,-5))
    @test all(BSE_SC.OneToIndex_to_Freq(10,10,5,5,0) .== (4,4))
    #shift tests
    @test all(BSE_SC.OneToIndex_to_Freq(1,1,5,5,1) .== (-5,-3))
    @test all(BSE_SC.OneToIndex_to_Freq(10,10,5,5,1) .== (4,2))
end

#TODO: reactivate, once github CI works
#= @testset "calc_χ₀" begin 
    gi_s0 = read_gImp("test_data/ED_s0.jld2")
    χ₀_s0 = calc_χ₀(gi_s0, 10.0, 20, 20, 0)
    gi_s1 = read_gImp("test_data/ED_s1.jld2")
    χ₀_s1 = calc_χ₀(gi_s1, 10.0, 20, 20, 1)
    @test χ₀_s0[0+21,0+21] ≈ -gi_s0[0]*gi_s0[0]*10.0
    @test χ₀_s1[0+21,0+21] ≈ -gi_s1[0]*gi_s1[0]*10.0
    @test χ₀_s0[10+21,5+21] ≈ -gi_s0[10]*gi_s0[5+10]*10.0
    @test χ₀_s1[10+21+2,5+21] ≈ -gi_s1[10]*gi_s1[5+10]*10.0
    @test χ₀_s0[10+21,-5+21] ≈ -gi_s0[10]*gi_s0[-5+10]*10.0
    @test χ₀_s1[10+21-2,-5+21] ≈ -gi_s1[10]*gi_s1[-5+10]*10.0
end
=#

@testset "n to i" begin
    @test BSE_SC.ωn_to_ωi(0) == 1
    @test BSE_SC.ωn_to_ωi(1) == 2
    @test BSE_SC.ωn_to_ωi(-1) == 2
    @test all(BSE_SC.OneToIndex_to_Freq(11,11,10,10,1) == (0,0))
    @test all(BSE_SC.OneToIndex_to_Freq(1,1,10,10,1) == (-10,-5))
    @test all(BSE_SC.OneToIndex_to_Freq(1,1,10,10,0) == (-10,-10))
end
