using Test
using BSE_SC

@testset "unit tests" begin
    @testset "IO" begin
        include("IO.jl")
    end

    @testset "IO" begin
        include("ladderDGA_core.jl")
    end

    @testset "chi asympt" begin
        include("build_chi_asympt.jl")
    end
end

@testset "functional tests" begin
    @testset "s0" begin
        include("test_s0.jl") 
    end
end
