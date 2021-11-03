using BSE_SC

NBorder = 10
gImp, χ₀, χDMFTsp_impr, χDMFTch_impr, χ_sp_asympt, χ_ch_asympt, χ_pp_asympt, U, β, shift = BSE_SC.setup("test_data/ED_s0.jld2", NBorder);
BSE_SC.improve_χ!(χDMFTsp_impr, χDMFTch_impr, χ₀, χ_sp_asympt, χ_ch_asympt, χ_pp_asympt, U, β, NBorder, shift);
