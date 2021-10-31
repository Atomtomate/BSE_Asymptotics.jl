@testset "types" begin
    include("types.jl")
end

@testset "chi" begin
    include("chi.jl")
end

@testset "chi asympt" begin
    include("build_chi_asympt.jl")
end
