struct BSE_SC_Helper
    χsp_asympt::Array{ComplexF64,1}
    χch_asympt::Array{ComplexF64,1}
    χpp_asympt::Array{ComplexF64,1}
    Fr::Array{ComplexF64,2}
    λr::Array{ComplexF64,1}
    Nν_shell::Int
    I_core::Array{CartesianIndex{2},1}
    I_corner::Array{CartesianIndex{2},1}
    I_t::Array{CartesianIndex{2},1}
    I_r::Array{CartesianIndex{2},1}
    I_all::Array{CartesianIndex{2},1}
    I_asympt::Array{CartesianIndex{2},1}
    ind1_list_corner::Array{Int,2}
    ind2_list_corner::Array{Int,2}

    function BSE_SC_Helper(χsp_asympt, χch_asympt, χpp_asympt, Nν_full, Nν_shell, n_iω, n_iν, shift)
        I_core, I_corner, I_t, I_r = shell_indices(Nν_full, Nν_shell)
        I_all = sort(union(I_core, I_corner, I_r, I_t))
        I_asympt = sort(union(I_corner, I_r, I_t))
        ind1_list_corner = Array{Int, 2}(undef, length(I_corner), 2*n_iω+1)
        ind2_list_corner = Array{Int, 2}(undef, length(I_corner), 2*n_iω+1)
        for ωi in 1:(2*n_iω+1)
            i1, i2 = aux_indices(I_corner, ωi, n_iω, n_iν, shift)
            ind1_list_corner[:,ωi] = i1
            ind2_list_corner[:,ωi] = i2
        end
        Fr = Array{ComplexF64,2}(undef, Nν_full, Nν_full)
        λr = Array{ComplexF64,1}(undef, Nν_full)

        new(χsp_asympt, χch_asympt, χpp_asympt, Fr, λr, Nν_shell, I_core, I_corner,
                            I_t, I_r, I_all, I_asympt, ind1_list_corner, ind2_list_corner)
    end
end

