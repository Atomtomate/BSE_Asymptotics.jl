"""
    shell_indices(Nν_full::Int, Nν_shell::Int)

Constructs 4 sets of indices for a square array of size `Nν_full` with a shell `Nν_shell` entries on each 
side in which asymptotic behaviour is assumed.
Returns `I_core`, `I_corner`, `I_top_bottom`, `I_left_right` with indices of the corresponding regions.
"""
function shell_indices(Nν_full::Int, Nν_shell::Int)
    ind_inner = (Nν_shell+1):(Nν_full-Nν_shell)
    ind_outer = union(1:(Nν_shell),(Nν_full-Nν_shell+1):Nν_full)
    I_top_bottom = sort([CartesianIndex(νi,νpi) for νi in ind_outer for νpi in ind_inner])
    I_left_right = sort([CartesianIndex(νi,νpi) for νi in ind_inner for νpi in ind_outer])
    I_corner = sort([CartesianIndex(νi,νpi) for νi in ind_outer for νpi in ind_outer])
    I_core = sort([CartesianIndex(νi,νpi) for νi in ind_inner for νpi in ind_inner])
    return I_core, I_corner, I_top_bottom, I_left_right
end

"""
    aux_indices(ind_lst::Vector{CartesianIndex{2}}, ωi::Int, n_iω::Int, n_iν::Int, shift::Int)

Constructs two index lists, one for ν-ν' and one for ν+ν'+ω, from input index list `ind_lst`.
These lists are used by [`improve_χ!`](@ref) interally.
"""
function aux_indices(ind_lst::Vector{CartesianIndex{2}}, ωi::Int, n_iω::Int, n_iν::Int, shift::Int)
    ind1_list = Vector{Int}(undef, length(ind_lst))    # ν' - ν
    ind2_list = Vector{Int}(undef, length(ind_lst))    # ν + ν' + ω
    for (i,ind) in enumerate(ind_lst)
        ind1_list[i] = ωn_to_ωi(ind[1] - ind[2])
        ωn, νn, = OneToIndex_to_Freq(ωi, ind[1], n_iω, n_iν, shift)
        ωn, νpn = OneToIndex_to_Freq(ωi, ind[2], n_iω, n_iν, shift)
        ind2_list[i] = ωn_to_ωi(νn + νpn + 1 + ωn)
    end
    return ind1_list, ind2_list
end

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

Generates helper struct, for the efficient computation of the susceptibility asymptotics,
using self consistency. See `BSE_Asym_Helper` for the helper for a direct version.
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

"""
    BSE_Asym_Helper

Helper struct for efficient solution of the self consistency. For an example see
also [`setup`](@ref).

Fields
-------------
- **`χsp_asympt`** : ω-asymptotic for the susceptibility in the spin-ph channel
- **`χch_asympt`** : ω-asymptotic for the susceptibility in the charge-ph channel
- **`χpp_asympt`** : ω-asymptotic for the susceptibility in the pp channel
- **`χ₀_shell_sum`** : Asymptotic sum for χ₀ shell
- **`Nν_shell`**   : Number of Fermionic frequencies used for asymptotic extension.
- **`I_core`**     : Indices for the core (non-asymptotic) region.
- **`I_asympt`**   : Indices for the asymptotic region (union of I_corner, I_t, I_r).
- **`ind1_list`**  : ν-ν' indices in `χsp_asympt` and `χch_asympt` for all `I_asympt`
- **`ind2_list`**  : ν+ν'+ω indices in `χpp_asympt` for all `I_asympt`
- **`shift`**      : `1` or `0` depending on wheter the ν-frequencies are shfited by `-ω/2`
"""
struct BSE_Asym_Helper
    χsp_asympt::Array{ComplexF64,1}
    χch_asympt::Array{ComplexF64,1}
    χpp_asympt::Array{ComplexF64,1}
    Nν_shell::Int
    I_core::Array{CartesianIndex{2},1}
    I_asympt::Array{CartesianIndex{2},1}
    ind1_list::Array{Int,1}
    ind2_list::OffsetArray{Int,2}
    shift::Int
    shell_sum_core::OffsetArray{ComplexF64,2}
    diag_asym_buffer::Array{ComplexF64,1}
"""
    BSE_Asym_Helper(χsp_asympt, χch_asympt, χpp_asympt, Nν_full, Nν_shell, n_iω, n_iν, shift)

Generates helper struct, for the efficient computation of the susceptibility asymptotics,
using self consistency. See `BSE_Asym_Helper` for the helper for a direct version.
`Nν_full`: number of Fermionic frequencies in ν direction for each Bosonic frequency.
`Nν_shell`: number of Fermionic frequencies in ν direction for each Bosonic frequency
`n_iω`: number of positive Bosonic frequencies.
`n_iν`: number of positive Fermionic frequencies.
`shift`: `1` or `0`, depending on whether or not the Fermionic frequencies are shifted by `ω/2`
"""
    function BSE_Asym_Helper(χsp_asympt, χch_asympt, χpp_asympt, Nν_shell, 
                             β, n_iω, n_iν, shift)
        Nν_full = 2*n_iν 
        I_core, I_corner, I_t, I_r = shell_indices(Nν_full, Nν_shell)
        I_all = sort(union(I_core, I_corner, I_r, I_t))
        I_asympt = sort(union(I_corner, I_r, I_t))
        ind2_list = OffsetArray(Array{Int, 2}(undef, length(I_asympt), 2*n_iω+1), 1:length(I_asympt), -n_iω:n_iω)
        i1, i2 = aux_indices(I_asympt, 1, n_iω, n_iν, shift)
        ind1_list = i1
        for ωi in 1:(2*n_iω+1)
            i1, i2 = aux_indices(I_asympt, ωi, n_iω, n_iν, shift)
            ind2_list[:,ωi-n_iω-1] = i2
        end
        shell_sum_core = χ₀_shell_sum_core(β, -n_iω:n_iω, n_iν - Nν_shell, shift)
        buffer = Array{ComplexF64,1}(undef, Nν_full)
        new(χsp_asympt, χch_asympt, χpp_asympt, Nν_shell, I_core, I_asympt,
            ind1_list, ind2_list, shift, shell_sum_core, buffer)
    end
end


"""
    χ₀_shell_sum_core(β::Float64, ω_ind_grid::AbstractVector{Int}, n_iν::Int, shift::Int)

Calculates the core region for use in [`χ₀_shell_sum`](@ref).
"""
function χ₀_shell_sum_core(β::Float64, ω_ind_grid::AbstractVector{Int}, n_iν::Int, shift::Int)
    res = OffsetArray(zeros(ComplexF64, length(ω_ind_grid), 4), ω_ind_grid, 1:4)
    iω_grid = iω_array(β, ω_ind_grid)
    for ωn in ω_ind_grid
        si = shift*trunc(Int,ωn/2)
        iν_ind_grid = (-n_iν+minimum(ω_ind_grid)-si):(n_iν+maximum(ω_ind_grid)-si-1)
        iν_grid = OffsetArray(iν_array(β, iν_ind_grid), iν_ind_grid)
        for νn in (-n_iν-si):(n_iν-si-1)
            res[ωn,1] += -1/(iν_grid[νn]^1 * iν_grid[νn+ωn]^1)
            res[ωn,2] +=  1/(iν_grid[νn]^2 * iν_grid[νn+ωn]^1) + 1/(iν_grid[νn]^1 * iν_grid[νn+ωn]^2)
            res[ωn,3] += -1/(iν_grid[νn]^3 * iν_grid[νn+ωn]^1) - 1/(iν_grid[νn]^1 * iν_grid[νn+ωn]^3)
            res[ωn,4] += -1/(iν_grid[νn]^2 * iν_grid[νn+ωn]^2)
        end
    end
    return res
end

"""
    χ₀_shell_sum(core::OffsetMatrix{ComplexF64}, ωn::Int, c2::Float64, c3::Float64)

Calculates the asymptotic sum of `∑ₙ,ₗ χ₀(iωₘ,iνₙ,iνₗ')` with `n,l` from `n_iν` to `∞`.
This is done by using the known first three tail coefficients of the Green's function in `iνₙ` and
expanding around `n → ∞`.
The core region is precalculated using [`χ₀_shell_sum_core`](@ref).
"""
function χ₀_shell_sum(core::OffsetArray{ComplexF64,2}, ωn::Int, β::Float64, c1::Float64, c2::Float64)
    res = (core[ωn,1] + c1*core[ωn,2] + c2*core[ωn,3] + c1^2*core[ωn,4])/β
    res += ωn == 0  ? -(β/4-c2*β^3/24-c1^2*β^3/48) : ((c1^2-c2)*β/(2*(2*ωn*π/β)^2))
    return res
end