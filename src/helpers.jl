abstract type BSE_Asym_Helpers end
abstract type BSE_SC_Helpers end

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
struct BSE_SC_Helper <: BSE_SC_Helpers
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
using self consistency. [`BSE_Asym_Helper`](@ref) is a helper for a direct asymptotics version.
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

Helper struct for efficient improvement of the result after the inversion in the Bethe Salpeter equation. 
For an example see also [`setup`](@ref).
TODO: update documentation `I_asympt` no longer in use. Optimize 

Fields
-------------
- **`χsp_asympt`** : ω-asymptotic for the susceptibility in the spin-ph channel
- **`χch_asympt`** : ω-asymptotic for the susceptibility in the charge-ph channel
- **`χpp_asympt`** : ω-asymptotic for the susceptibility in the pp channel
- **`Nν_shell`**   : Number of Fermionic frequencies used for asymptotic extension.
- **`I_core`**     : Indices for the core (non-asymptotic) region.
- **`block_i`**    : Sorted index into `F_asympt`.
- **`block_slices`** : Blocks of indices for the remaining the index lists (this improves performance by collecting all summands for one F_asymp index into one block).
- **`ind1_list`**  : index into χ₀ for all `I_asympt
- **`ind2_list`**  : ν-ν' indices in `χsp_asympt` and `χch_asympt` for all `I_asympt`
- **`ind3_list`**  : ν+ν'+ω indices in `χpp_asympt` for all `I_asympt`
- **`shift`**      : `1` or `0` depending on wheter the ν-frequencies are shfited by `-ω/2`
"""
struct BSE_Asym_Helper <: BSE_Asym_Helpers
    χsp_asympt::Array{ComplexF64,1}
    χch_asympt::Array{ComplexF64,1}
    χpp_asympt::Array{ComplexF64,1}
    Nν_shell::Int
    I_core::Array{CartesianIndex{2},1}
    block_i::Vector{Int}
    block_slices::Vector{UnitRange{Int64}}
    ind1_list::Array{Int,1}
    ind2_list::Array{Int,1}
    ind3_list::OffsetArray{Int,2}
    shift::Int
    diag_asym_buffer::Array{ComplexF64,1}
"""
    BSE_Asym_Helper(χsp_asympt, χch_asympt, χpp_asympt, Nν_full, Nν_shell, n_iω, n_iν, shift)

Generates helper struct, for the efficient computation of the susceptibility asymptotics,
using direct asymptotics. See `BSE_Asym_Helper` for the helper for a direct version.
`Nν_full`: number of Fermionic frequencies in ν direction for each Bosonic frequency.
`Nν_shell`: number of Fermionic frequencies in ν direction for each Bosonic frequency
`n_iω`: number of positive Bosonic frequencies.
`n_iν`: number of positive Fermionic frequencies.
`shift`: `1` or `0`, depending on whether or not the Fermionic frequencies are shifted by `ω/2`
"""
    function BSE_Asym_Helper(χsp_asympt::Vector{ComplexF64}, 
                             χch_asympt::Vector{ComplexF64}, 
                             χpp_asympt::Vector{ComplexF64}, 
                             Nν_shell::Int, U::Float64, β::Float64, n_iω::Int, n_iν::Int, shift::Int)
        n_iν_f = n_iν + Nν_shell
        Nν_full = 2*n_iν_f 
        I_core, I_corner, I_t, I_r = shell_indices(Nν_full, Nν_shell)
        I_asympt = sort(union(I_corner, I_r, I_t))
        ind2_list = OffsetArray(Array{Int, 2}(undef, length(I_asympt), 2*n_iω+1), 1:length(I_asympt), -n_iω:n_iω)
        i1l, i2l = aux_indices(I_asympt, 1, n_iω, n_iν_f, shift)
        for ωi in 1:(2*n_iω+1)
            i1l, i2l = aux_indices(I_asympt, ωi, n_iω, n_iν_f, shift)
            ind2_list[:,ωi-n_iω-1] = i2l
        end
        block_i, block_slices, ind1_list, ii = optimize_indices(I_asympt)
        buffer = Array{ComplexF64,1}(undef, Nν_full)
        new(χsp_asympt, χch_asympt, χpp_asympt, Nν_shell, I_core, block_i, block_slices, ind1_list,
            i1l[ii], ind2_list[ii,:], shift, buffer)
    end
end

"""
    optimize_indices(I_asympt)

Returns blocks of consecutive indices for the first index in `I_asympt` and corresponding slices for the indexing of the remaining indices.
`block_indices` contains the slices of indices used in the remaining three index lists (`ind1_list` and `ind2_list` in `BSE_Asym_Helper`).
`ii_2[ii]` is the second index of the `I_asympt` list.
"""
function optimize_indices(I_asympt)
    ii_1::Vector{Int} = map(x->x[1], I_asympt)
    ii_2::Vector{Int} = map(x->x[2], I_asympt)

    ii = sortperm(ii_1)
    tmp_i1_l = ii_1[ii]
    block_slices::Vector{UnitRange{Int64}} = []
    block_i::Vector{Int} = []
    ci = 1
    ci_s = 1
    for i in eachindex(tmp_i1_l)
        if tmp_i1_l[i] > ci
            push!(block_slices, ci_s:(i-1))
            push!(block_i, ci)
            ci = tmp_i1_l[i]
            ci_s = i
        end
    end
    push!(block_slices, ci_s:length(tmp_i1_l))
    push!(block_i, tmp_i1_l[end])
    block_i, block_slices, ii_2[ii], ii
end

"""
    BSE_Asym_Helper_Approx1

Helper struct for efficient improvement of the result after the inversion in the Bethe Salpeter equation. 
This approximates the results by neglecting the partice-particle channel, greatly imrpoving performance.
In rare cases where the particle-particle susceptibility contributes significantly use [`BSE_Asym_Helper`](@ref) instead.

Fields
-------------
- **`Nν_shell`**   : Number of Fermionic frequencies used for asymptotic extension.
- **`Fsp_kernel`** : Approximate diagonal contribution from physical susceptibilities for the spin channel
- **`Fch_kernel`** : Approximate diagonal contribution from physical susceptibilities for the charge channel
"""
struct BSE_Asym_Helper_Approx1 <: BSE_Asym_Helpers
    Nν_shell::Int
    Fsp_kernel::Matrix{ComplexF64}
    Fch_kernel::Matrix{ComplexF64}
    diag_asym_buffer::Array{ComplexF64,1}
"""
    BSE_Asym_Helper_Approx1(χsp_asympt, χch_asympt, χpp_asympt, Nν_full, Nν_shell, n_iω, n_iν, shift)

Generates helper struct, for the efficient computation of the susceptibility asymptotics,
using direct asymptotics. 
Note that this approximation disregards the `χpp_asympt` input!
`Nν_full`: number of Fermionic frequencies in ν direction for each Bosonic frequency.
`Nν_shell`: number of Fermionic frequencies in ν direction for each Bosonic frequency
`n_iω`: number of positive Bosonic frequencies.
`n_iν`: number of positive Fermionic frequencies.
`shift`: `1` or `0`, depending on whether or not the Fermionic frequencies are shifted by `ω/2`
"""
    function BSE_Asym_Helper_Approx1(χsp_asympt, χch_asympt, χpp_asympt, Nν_shell, U, β, n_iω, n_iν, shift)
        n_iν_f = n_iν + Nν_shell
        Nν_full = 2*n_iν_f 
        I_core, I_corner, I_t, I_r = shell_indices(Nν_full, Nν_shell)
        I_all = sort(union(I_core, I_corner, I_r, I_t))
        I_asympt = sort(union(I_corner, I_r, I_t))
        ind2_list = OffsetArray(Array{Int, 2}(undef, length(I_asympt), 2*n_iω+1), 1:length(I_asympt), -n_iω:n_iω)
        i1l, i2l = aux_indices(I_asympt, 1, n_iω, n_iν_f, shift)
        Fsp_kernel = Matrix{ComplexF64}(undef, 2*n_iν, 2*n_iν)
        Fch_kernel = Matrix{ComplexF64}(undef, 2*n_iν, 2*n_iν)

        for ωi in 1:(2*n_iω+1)
            i1l, i2l = aux_indices(I_asympt, ωi, n_iω, n_iν_f, shift)
            for i in 1:length(i1l)
                i1 = I_asympt[i]
                i2 = i1l[i]
            end
        end

        buffer = Array{ComplexF64,1}(undef, Nν_full)
        #println(typeof(χsp_asympt), typeof(χch_asympt), typeof(χpp_asympt), typeof(χsp_asym_b), typeof(χch_asym_b), typeof(Nν_shell), typeof(I_core), typeof(I_asympt),
        new(Nν_shell, Fsp_kernel, Fch_kernel, buffer)
    end
end

"""
    BSE_Asym_Helper_Approx2

Helper struct for efficient improvement of the result after the inversion in the Bethe Salpeter equation. 
This approximates the results by neglecting all offdiagonal terms greatly imrpoving performance and removing dependence on the physical susceptibilities.

Fields
-------------
- **`Nν_shell`**   : Number of Fermionic frequencies used for asymptotic extension.
"""
struct BSE_Asym_Helper_Approx2 <: BSE_Asym_Helpers
    Nν_shell::Int
"""
    BSE_Asym_Helper_Approx2(Nν_full, Nν_shell, n_iω, n_iν, shift)

Generates helper struct, for the efficient computation of the susceptibility asymptotics,
using direct asymptotics. 
Note that this approximation disregards the `χpp_asympt` input!
`Nν_shell`: number of Fermionic frequencies in ν direction for each Bosonic frequency
"""
    function BSE_Asym_Helper_Approx2(Nν_shell)
        new(Nν_shell)
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
        iν_ind_grid = (-n_iν-abs(minimum(ω_ind_grid))-si):(n_iν+abs(maximum(ω_ind_grid))-si-1)
        iν_grid = OffsetArray(iν_array(β, iν_ind_grid), iν_ind_grid)
        for νn in (-n_iν-si):(n_iν-si-1)
            res[ωn,1] += 1/(iν_grid[νn]^1 * iν_grid[νn+ωn]^1)
            res[ωn,2] += 1/(iν_grid[νn]^2 * iν_grid[νn+ωn]^1) + 1/(iν_grid[νn]^1 * iν_grid[νn+ωn]^2)
            res[ωn,3] += 1/(iν_grid[νn]^2 * iν_grid[νn+ωn]^2)
            res[ωn,4] += 1/(iν_grid[νn]^3 * iν_grid[νn+ωn]^1) + 1/(iν_grid[νn]^1 * iν_grid[νn+ωn]^3)
        end
    end
    return res
end

"""
    χ₀_shell_sum(core::OffsetMatrix{ComplexF64}, ωn::Union{Int,AbstractVector{Int}}, c2::Float64, c3::Float64)

Calculates the asymptotic sum of `∑ₙ,ₗ χ₀(iωₘ,iνₙ,iνₗ')` with `n,l` from `n_iν` to `∞`.
This is done by using the known first three tail coefficients of the Green's function in `iνₙ` and
expanding around `n → ∞`.
`ωn` can either be a single Matsubara index or an array of indices. 
The core region is precalculated using [`χ₀_shell_sum_core`](@ref).
"""
function χ₀_shell_sum(core::OffsetArray{ComplexF64,2}, ωn_in::Int, β::Float64, c1::Float64, c2::Float64, c3::Float64)::ComplexF64
    @inbounds res = (core[ωn_in,1] + c1*core[ωn_in,2] + c2*core[ωn_in,3] + c3*core[ωn_in,4])/β
    res += ωn_in == 0 ? -(-β/4+c2*β^3/48+c3*β^3/24) : ((c3-c2)*β/(2*(2*ωn_in*π/β)^2))
    return res
end

function χ₀_shell_sum(core::OffsetArray{ComplexF64,2}, ωn_in::AbstractVector{Int}, β::Float64, c1::Float64, c2::Float64, c3::Float64)::ComplexF64
    res = Array{ComplexF64,1}(undef, length(ωn_in))
    for (ωi,ωn) in enumerate(ωn_in)
        @inbounds res[ωi] = (core[ωn,1] + c1*core[ωn,2] + c2*core[ωn,3] + c3*core[ωn,4])/β
        res += ωn == 0 ? -(-β/4+c2*β^3/48+c3*β^3/24) : ((c3-c2)*β/(2*(2*ωn*π/β)^2))
    end
    return res
end
