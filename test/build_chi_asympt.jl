@testset "shell indices" begin 
    Nν_full = 7
    Nν_shell = 2
    test_arr = reshape(1:(Nν_full*Nν_full), (Nν_full,Nν_full))
    I_core, I_corner, I_top, I_side = BSE_SC.shell_indices(Nν_full, Nν_shell)
    @test all(test_arr[I_corner] .== [1,8,36,43,2,9,37,44,6,13,41,48,7,14,42,49])
    @test all(test_arr[I_top] .== [15,22,29,16,23,30,20,27,34,21,28,35])
    @test all(test_arr[I_side] .== [3,10,38,45,4,11,39,46,5,12,40,47])
    @test all(test_arr[I_core] .== [17,24,31,18,25,32,19,26,33])
end

@testset "auxiliary indices" begin
    Nν_full = 7
    Nν_shell = 2
    shift = 0
    test_arr = reshape(1:(Nν_full*Nν_full), (Nν_full,Nν_full))
    I_core, I_corner, I_top, I_side = BSE_SC.shell_indices(Nν_full, Nν_shell)
    ωi = 3
    n_iω = 5
    n_iν = trunc(Int, Nν_full/2)
    ind1_c, ind2_c = BSE_SC.aux_indices(I_corner, ωi, n_iω, n_iν, shift)
end