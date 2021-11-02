gi = read_gImp("test_data/ED_out_b10.0_u1.0_20_20_s0.jld2")
@test gi.offsets[1] == -83
@test gi[0] ≈ 2.0e-16 - 1.3645385209655727im
@test gi[-1] ≈ conj(gi[0])
@test gi[-10] ≈ conj(gi[9])
