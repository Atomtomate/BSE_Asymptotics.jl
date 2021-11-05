using LinearAlgebra
using JLD2
using OffsetArrays
using PaddedViews

testf = joinpath(@__DIR__,"test_data/ED_out_large.jld2")
testf_s1 = joinpath(@__DIR__,"test_data/ED_s1.jld2")
testf_s0 = joinpath(@__DIR__,"test_data/ED_s0.jld2")
include(joinpath(@__DIR__,"../src/LapackWrapper.jl"))
include(joinpath(@__DIR__,"../src/ladderDGA_core.jl"))
include(joinpath(@__DIR__,"../src/IO.jl"))
include(joinpath(@__DIR__,"../src/build_chi_asympt.jl"))


NBorder = 20


gImp_s0, χ₀_s0, χDMFTsp_impr_s0, χDMFTch_impr_s0, χ_sp_asympt_s0, χ_ch_asympt_s0, χ_pp_asympt_s0, U_s0, β_s0, shift_s0 = 
        setup(testf_s0, NBorder);
χDMFTsp_s0 = deepcopy(χDMFTsp_impr_s0)
χDMFTch_s0 = deepcopy(χDMFTch_impr_s0)

improve_χ!(χDMFTsp_impr_s0, χDMFTch_impr_s0, χ₀_s0, χ_sp_asympt_s0, χ_ch_asympt_s0, χ_pp_asympt_s0, U_s0, β_s0, NBorder, shift_s0);

gImp_s1, χ₀_s1, χDMFTsp_impr_s1, χDMFTch_impr_s1, χ_sp_asympt_s1, χ_ch_asympt_s1, χ_pp_asympt_s1, U_s1, β_s1, shift_s1 = 
        setup(testf_s1, NBorder);
χDMFTsp_s1 = deepcopy(χDMFTsp_impr_s1)
χDMFTch_s1 = deepcopy(χDMFTch_impr_s1)

improve_χ!(χDMFTsp_impr_s1, χDMFTch_impr_s1, χ₀_s1, χ_sp_asympt_s1, χ_ch_asympt_s1, χ_pp_asympt_s1, U_s1, β_s1, NBorder, shift_s1);