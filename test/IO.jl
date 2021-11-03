gi = read_gImp("test_data/ED_s0.jld2")
@test gi.offsets[1] == -83
@test gi[0] ≈ 2.0e-16 - 1.3645385209655727im
@test gi[-1] ≈ conj(gi[0])
@test gi[-10] ≈ conj(gi[9])
