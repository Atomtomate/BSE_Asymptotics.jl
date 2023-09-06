@testset "BSE_Asymptotics" begin
    @testset "shell indices" begin 
        Nν_full = 7
        Nν_shell = 2
        test_arr = reshape(1:(Nν_full*Nν_full), (Nν_full,Nν_full))
        I_core, I_corner, I_top, I_side = BSE_Asymptotics.shell_indices(Nν_full, Nν_shell)
        @test all(test_arr[I_corner] .== sort([1,8,36,43,2,9,37,44,6,13,41,48,7,14,42,49]))
        @test all(test_arr[I_top] .== sort([15,22,29,16,23,30,20,27,34,21,28,35]))
        @test all(test_arr[I_side] .== sort([3,10,38,45,4,11,39,46,5,12,40,47]))
        @test all(test_arr[I_core] .== sort([17,24,31,18,25,32,19,26,33]))
    end

    @testset "auxiliary indices" begin
        Nν_full = 8
        Nν_shell = 2
        shift = 0
        I_core, I_corner, I_top, I_side = BSE_Asymptotics.shell_indices(Nν_full, Nν_shell)
        n_iω = 5
        n_iν = trunc(Int, Nν_full/2)
        ind1_check = [1,2,     7,8,
                      2,1,     6,7,

                      7,6,     1,2,
                      8,7,     2,1]
        ind2_check = [8,7,     2,1,
                      7,6,     1,2,

                      2,1,     6,7,
                      1,2,     7,8]
        ind1_s0, ind2_s0 = BSE_Asymptotics.aux_indices(I_corner, 1, n_iω, n_iν, shift)
        @test all(ind1_s0 .== ind1_check)
        ind1_s0, ind2_s0 = BSE_Asymptotics.aux_indices(I_corner, 3, n_iω, n_iν, shift)
        @test all(ind1_s0 .== ind1_check)
        ind1_s0, ind2_s0 = BSE_Asymptotics.aux_indices(I_corner, 4, n_iω, n_iν, shift)
        @test all(ind1_s0 .== ind1_check)
    end
end

@testset "BSE_Asym" begin
    h = BSE_Asym_Helper(1:20,100:120,1000:1020,2,1.1,11.2,3,3,1)
    h_a1 = BSE_Asym_Helper_Approx1(1:20,100:120,1000:1020,2,1.1,11.2,3,3,1)
    h_a2 = BSE_Asym_Helper_Approx2(3)
end
