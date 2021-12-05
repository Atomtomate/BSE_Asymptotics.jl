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

#TODO: index struct
"""
    improve_χ!(type::Symbol, ωi::Int, χr, χ₀, U, β, h::BSE_SC_Helper; 
               Nit=200, atol=1e-9)

Improves asymptotics of `χr`, given a channel `type = :sp` or `:ch`, using Eq. 12a 12b 
from DOI: 10.1103/PhysRevB.97.235.140
`ωi` specifies the index of the Matsubara frequency to use. `χ₀` is the bubble term,
`U` and `β` Hubbard-U and temperature.
`h` is a helper struct, see [`BSE_SC_Helper`](@ref).
Additionally one can specify convergence parameters:
`Nit`:  maximum number of iterations
`atol`: minimum change between total sum over `χ` between iterations. 

"""
function improve_χ!(type::Symbol, ωi::Int, χr::AbstractArray{ComplexF64,2}, χ₀::AbstractArray{ComplexF64,1}, 
                U::Float64, β::Float64, h::BSE_SC_Helper; Nit=200, atol=1e-9)
    f = if type == :sp
        update_Fsp!
    elseif type == :ch
        update_Fch!
    else
        error("Unkown channel. Only sp/ch implemented")
    end

    fill!(h.Fr, 0.0)
    fill!(h.λr, 0.0)
    for i in h.I_core
        δ_ννp = Float64(i[1] == i[2])
        h.Fr[i] = - β^2 * (χr[i] - δ_ννp*χ₀[i[1]])/(χ₀[i[1]]*χ₀[i[2]])
    end
    χr_old = 0.0
    converged = false
    i = 0
    while !converged && (i < Nit)
        χr_n = update_χ!(h.λr, χr, h.Fr, χ₀, β, h.I_asympt)
        f(χr_n, U, ωi, h)
        if (abs(χr_old - χr_n) < atol)
            converged = true
        else
            i += 1
            χr_old = χr_n
        end
    end
    return i
end

function update_Fsp!(χ::ComplexF64, U::Float64, ωi::Int, h::BSE_SC_Helper)
    i1_l = h.ind1_list
    i2_l = view(h.ind2_list, :, ωi)
    for i in 1:length(i1_l)
        i1 = h.I_asympt[i]
        i2 = i1_l[i]
        i3 = i2_l[i]
        h.Fr[i1] = -U + (U^2)*χ + U*h.λr[i1[1]]  + U*h.λr[i1[2]] 
            + (U^2/2)*h.χch_asympt[i2] - (U^2/2)*h.χsp_asympt[i2] + (U^2)*h.χpp_asympt[i3] 
    end
end

function update_Fch!(χ::ComplexF64, U::Float64, ωi::Int, h::BSE_SC_Helper)
    i1_l = h.ind1_list
    i2_l = view(h.ind2_list, :, ωi)
    for i in 1:length(i1_l)
        i1 = h.I_asympt[i]
        i2 = i1_l[i]
        i3 = i2_l[i]
        h.Fr[i1] = U + (U^2)*χ  - U*h.λr[i1[1]] - U*h.λr[i1[2]]
            + (U^2/2)*h.χch_asympt[i2] + 3*(U^2/2)*h.χsp_asympt[i2] - (U^2)*h.χpp_asympt[i3]
    end
end

function update_χ!(λ::AbstractArray{ComplexF64,1}, χ::AbstractArray{ComplexF64,2}, 
        F::AbstractArray{ComplexF64,2}, χ₀::AbstractArray{ComplexF64,1}, 
        β::Float64, indices::AbstractArray{CartesianIndex{2},1})
    for i in indices
        x[i] = i[1] == i[2] ? χ₀[i[1]] : 0
        χ[i] -= χ₀[i[1]]*F[i]*χ₀[i[2]]/(β^2)
    end
    for νk in 1:size(χ,2)
        λ[νk] = (sum(view(χ,:,νk)) / (-χ₀[νk])) + 1
    end
    return sum(χ)/(β^2)
end

function χ₀sum(iνₙ)
    
end

function improve_χλ_direct(χsp::AbstractArray{ComplexF64,2}, χch::AbstractArray{ComplexF64,2}, χ₀::AbstractArray{ComplexF64,1}, U::Float64, β::Float64, bs, h::BSE_SC_Helper)
    i1_l = h.ind1_list
    i2_l = view(h.ind2_list, :, ωi)
    χ₀_core = view(χ₀,(h.Nν_shell+1):(length(χ₀)-h.Nν_shell))
    λsp_core = -sum(χsp,dims=[2])[:,1] ./ χ₀_core .+ 1
    λch_core =  sum(χch,dims=[2])[:,1] ./ χ₀_core .- 1
    χsp_core = sum(χsp) /β^2
    χch_core = sum(χch) /β^2
    λasym_ch = zeros(eltype(h.χch_asympt), length(h.λr))
    λasym_sp = zeros(eltype(h.χch_asympt), length(h.λr))
    #TODO this is inefficient. only loop over necessary indices, only one direction needs to be extended
    for i in 1:length(i1_l)
        i1 = h.I_asympt[i]
        i2 = i1_l[i]
        i3 = i2_l[i]
        λasym_ch[i1[1]] += ((U^2/2)*h.χch_asympt[i2] + 3*(U^2/2)*h.χsp_asympt[i2] - (U^2)*h.χpp_asympt[i3])*(-χ₀[i1[2]])/β^2
        λasym_sp[i1[1]] += ((U^2/2)*h.χch_asympt[i2] -   (U^2/2)*h.χsp_asympt[i2] + (U^2)*h.χpp_asympt[i3])*(-χ₀[i1[2]])/β^2
    end
    λasym_sp_core = view(λasym_sp,(h.Nν_shell+1):(length(χ₀)-h.Nν_shell))
    λasym_ch_core = view(λasym_ch,(h.Nν_shell+1):(length(χ₀)-h.Nν_shell))
    λsp = (λsp_core .+ λasym_sp_core .+ U*bs)/(1+U*bs)
    λch = (λch_core .+ λasym_ch_core .+ U*bs)/(1-U*bs)
    λsp_s = -sum((λsp .- 1) .* χ₀_core)/β^2
    λch_s = -sum((λch .+ 1) .* χ₀_core)/β^2
    λasym_sp_s = -sum(λasym_sp_core .* χ₀_core)/β^2
    λasym_ch_s = -sum(λasym_ch_core .* χ₀_core)/β^2
    χsp = (χsp_core - bs*(1+U*(2*λsp_s-bs)) - U^2 * λasym_sp_s)/(1-U^2 * bs^2)
    χch = (χch_core - bs*(1+U*(2*λch_s+bs)) - U^2 * λasym_ch_s)/(1-U^2 * bs^2)
    return χsp, χch, λsp, λch
end

function χ_direct(χ::AbstractArray{ComplexF64,2}, χ₀::AbstractArray{ComplexF64,1}, U::Float64, β::Float64, h::BSE_SC_Helper)
    i1_l = h.ind1_list
    i2_l = view(h.ind2_list, :, ωi)
    χconst_ch = zeros(eltype(h.χch_asympt), length(i1_l))
    χconst_sp = zeros(eltype(h.χch_asympt), length(i1_l))
    χ₀_asym = χ₀[union(1:10,111:120)]
    χ₀_core = χ₀[11:110]
    f1 = 1 / (1 - (U^2 / β^4) * sum(χ₀ .* χ₀))
    χ_ch = sum(χ₀)/β^2 - sum(χ)/β^2
    χ_sp = sum(χ₀)/β^2 - sum(χ)/β^2#sum(χ₀_core)/(β) + U*sum(χ₀_asym .* χ₀_asym)/(β^2) - sum(χ)/β^2
    for i in 1:length(i1_l)
        i1 = h.I_asympt[i]
        i2 = i1_l[i]
        i3 = i2_l[i]
        χconst_ch[i] = χ₀[i1[1]]*((U^2/2)*h.χch_asympt[i2] + 3*(U^2/2)*h.χsp_asympt[i2] - (U^2)*h.χpp_asympt[i3])*χ₀[i1[2]]
        χconst_sp[i] = χ₀[i1[1]]*((U^2/2)*h.χch_asympt[i2] -   (U^2/2)*h.χsp_asympt[i2] + (U^2)*h.χpp_asympt[i3])*χ₀[i1[2]]
    end
    χ_ch -= sum(χconst_ch)/β^4
    χ_sp -= sum(χconst_sp)/β^4
    return χ_sp,χ_ch,f1*χ_sp, f1*χ_ch
end
