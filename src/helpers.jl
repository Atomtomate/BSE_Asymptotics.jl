"""
    BSE_SC_Helper

Helper struct for efficient solution of the self consistency. For an example see
also [`setup`](@ref).

Fields
-------------
- **`χsp_asympt`** : ω-asymptotic for the susceptibility in the spin-ph channel
- **`χch_asympt`** : ω-asymptotic for the susceptibility in the charge-ph channel
- **`χpp_asympt`** : ω-asymptotic for the susceptibility in the pp channel
- **`Fr`**         : Temporary storage for the full vertex. Used internally to avoid repeated mallocs.
- **`λr`**         : Temporary storage for the fermion-boson-vertex. Used internally to avoid repeated mallocs.
- **`Nν_shell`**   : Number of Fermionic frequencies used for asymptotic extension.
- **`I_core`**     : Indices for the core (non-asymptotic) region.
- **`I_corner`**   : Indices for the corner (ν AND ν' outside core) region.
- **`I_t`**        : Indices for the core (ν outside core) region.
- **`I_r`**        : Indices for the core (ν' outside core) region.
- **`I_asympt`**   : Indices for the asymptotic region (union of I_corner, I_t, I_r).
- **`ind1_list`**  : ν-ν' indices in `χsp_asympt` and `χch_asympt` for all `I_asympt`
- **`ind2_list`**  : ν+ν'+ω indices in `χpp_asympt` for all `I_asympt`
- **`shift`**      : `1` or `0` depending on wheter the ν-frequencies are shfited by `-ω/2`
"""
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
    I_asympt::Array{CartesianIndex{2},1}
    ind1_list::Array{Int,1}
    ind2_list::Array{Int,2}
    shift::Int

"""
    BSE_SC_Helper(χsp_asympt, χch_asympt, χpp_asympt, Nν_full, Nν_shell, n_iω, n_iν, shift)

Generates helper struct, for the efficient computation of the susceptibility asymptotics.
`Nν_full`: number of Fermionic frequencies in ν direction for each Bosonic frequency.
`Nν_shell`: number of Fermionic frequencies in ν direction for each Bosonic frequency
`n_iω`: number of positive Bosonic frequencies.
`n_iν`: number of positive Fermionic frequencies.
`shift`: `1` or `0`, depending on whether or not the Fermionic frequencies are shifted by `ω/2`
"""
    function BSE_SC_Helper(χsp_asympt, χch_asympt, χpp_asympt, Nν_full, Nν_shell, n_iω, n_iν, shift)
        #TODO: Nν_full = 2*Nν_shell+2*n_iν 
        I_core, I_corner, I_t, I_r = shell_indices(Nν_full, Nν_shell)
        I_all = sort(union(I_core, I_corner, I_r, I_t))
        I_asympt = sort(union(I_corner, I_r, I_t))
        ind2_list = Array{Int, 2}(undef, length(I_asympt), 2*n_iω+1)
        i1, i2 = aux_indices(I_asympt, 1, n_iω, n_iν, shift)
        ind1_list = i1
        for ωi in 1:(2*n_iω+1)
            i1, i2 = aux_indices(I_asympt, ωi, n_iω, n_iν, shift)
            ind2_list[:,ωi] = i2
        end
        Fr = Array{ComplexF64,2}(undef, Nν_full, Nν_full)
        λr = Array{ComplexF64,1}(undef, Nν_full)

        new(χsp_asympt, χch_asympt, χpp_asympt, Fr, λr, Nν_shell, I_core, I_corner,
                            I_t, I_r, I_asympt, ind1_list, ind2_list, shift)
    end
end

