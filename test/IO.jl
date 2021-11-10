#TODO: reactivate, once github CI works
#=gi = read_gImp("test_data/ED_s0.jld2")
@test gi[0] ≈ -3.8e-15 - 1.5455997882742398im
@test gi[-1] ≈ conj(gi[0])
@test gi[-10] ≈ conj(gi[9])
=#
