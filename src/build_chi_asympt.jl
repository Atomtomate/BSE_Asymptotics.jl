#TODO: index struct
"""
    improve_χ!(type::Symbol, ωi::Int, χr, χ₀, U, β, h::BSE_SC_Helper; 
               Nit=200, atol=1e-9)

Improves asymptotics of `χr`, given a channel `type = :sp` or `:ch`, using Eq. 12a 12b 
from DOI: 10.1103/PhysRevB.97.235.140
`ωn` specifies the index of the Matsubara frequency to use. `χ₀` is the bubble term,
`U` and `β` Hubbard-U and temperature.
`h` is a helper struct, see [`BSE_SC_Helper`](@ref).
Additionally one can specify convergence parameters:
`Nit`:  maximum number of iterations
`atol`: minimum change between total sum over `χ` between iterations. 

"""
function improve_χ!(type::Symbol, ωn::Int, χr::AbstractArray{ComplexF64,2}, χ₀::AbstractArray{ComplexF64,1}, 
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
        f(χr_n, U, ωn, h)
        if (abs(χr_old - χr_n) < atol)
            converged = true
        else
            i += 1
            χr_old = χr_n
        end
    end
    return i
end

function update_Fsp!(χ::ComplexF64, U::Float64, ωn::Int, h::BSE_SC_Helper)
    i1_l = h.ind1_list
    i2_l = view(h.ind2_list, :, ωn)
    for i in 1:length(i1_l)
        i1 = h.I_asympt[i]
        i2 = i1_l[i]
        i3 = i2_l[i]
        h.Fr[i1] = -U + (U^2)*χ + U*h.λr[i1[1]]  + U*h.λr[i1[2]] 
            + (U^2/2)*h.χch_asympt[i2] - (U^2/2)*h.χsp_asympt[i2] + (U^2)*h.χpp_asympt[i3] 
    end
end

function update_Fch!(χ::ComplexF64, U::Float64, ωn::Int, h::BSE_SC_Helper)
    i1_l = h.ind1_list
    i2_l = view(h.ind2_list, :, ωn)
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

function F_diag!(type::Symbol, ωn::Int, U::Float64, β::Float64, χ₀::AbstractArray{ComplexF64,1}, h::BSE_Asym_Helper)
    i1_l = h.ind1_list
    i2_l = view(h.ind2_list, :, ωn)
    fill!(h.diag_asym_buffer, 0)
    if type == :sp
        for i in 1:length(i1_l)
            i1 = h.I_asympt[i]
            i2 = i1_l[i]
            i3 = i2_l[i]
            h.diag_asym_buffer[i1[1]] += ((U^2/2)*h.χch_asympt[i2] - (U^2/2)*h.χsp_asympt[i2] + (U^2)*h.χpp_asympt[i3])*(-χ₀[i1[2]])/β^2
        end
    elseif type == :ch
        for i in 1:length(i1_l)
            i1 = h.I_asympt[i]
            i2 = i1_l[i]
            i3 = i2_l[i]
            h.diag_asym_buffer[i1[1]] += ((U^2/2)*h.χch_asympt[i2] + 3*(U^2/2)*h.χsp_asympt[i2] - (U^2)*h.χpp_asympt[i3])*(-χ₀[i1[2]])/β^2
        end
    else
        error("Unrecognized type $(type) for F_diag! Expected sp/ch")
    end
end

F_diag!(type, ωn, U, β, χ₀, h::BSE_Asym_Helper_Approx2) = nothing


"""
    calc_χλ(type::Symbol, ωn::Int, χ::AbstractArray{ComplexF64,2}, χ₀::AbstractArray{ComplexF64,1}, U::Float64, β::Float64, χ₀_asym::Float64, h::BSE_Asym_Helper)

Calculates the physical susceptibility `χ` and triangular vertex `λ` in a given channel `type=:sp` or `type=:ch` using knowledge about the asymptotics of the full vertex and tails of the Green's function.
`χ₀_asym` is the `ω` dependent asymptotic tail of `χ₀` and can be calculated with  [`χ₀_shell_sum`](@ref).

TODO: refactor code duplications
"""
function calc_χλ_impr(type::Symbol, ωn::Int, χ::AbstractArray{ComplexF64,2}, χ₀::AbstractArray{ComplexF64,1}, 
                 U::Float64, β::Float64, χ₀_asym::ComplexF64, h::HT) where  HT <: BSE_Asym_Helpers
    λ = Array{eltype(χ),1}(undef, size(χ, 1))
    calc_χλ_impr!(λ, type, ωn, χ, χ₀, U, β, χ₀_asym, h)
    return χ_out, λ
end

function calc_χλ_impr!(λ::Array{ComplexF64,1}, type::Symbol, ωn::Int, χ::AbstractArray{ComplexF64,2}, χ₀::AbstractArray{ComplexF64,1}, 
                 U::Float64, β::Float64, χ₀_asym::ComplexF64, h::BSE_Asym_Helper)
    s = type == :ch ? -1 : 1
    ind_core = (h.Nν_shell+1):(size(χ₀,1)-h.Nν_shell)
    χ₀_core = view(χ₀,ind_core)
    λ[:] = -s*sum(χ,dims=[2])[:,1] ./ χ₀_core .+ s
    χ_core = sum(χ)/β^2
    F_diag!(type, ωn, U, β, χ₀, h)
    λ[:] = (λ .- s*view(h.diag_asym_buffer, ind_core) .- U*χ₀_asym)/(1-s*U*χ₀_asym)
    λ_s = -sum((U .* λ .- s*U) .* χ₀_core)/β^2
    diag_asym_s = -sum(h.diag_asym_buffer .* χ₀)/β^2
    χ_out = (χ_core + χ₀_asym*(1+2*λ_s+s*U*χ₀_asym) - diag_asym_s)/(1-U^2 * χ₀_asym^2)
    return χ_out
end

function calc_χλ_impr!(λ::Array{ComplexF64,1}, type::Symbol, ωn::Int, χ::AbstractArray{ComplexF64,2}, χ₀::AbstractArray{ComplexF64,1}, 
                 U::Float64, β::Float64, χ₀_asym::ComplexF64, h::BSE_Asym_Helper_Approx2)
    s = type == :ch ? -1 : 1
    ind_core = (h.Nν_shell+1):(size(χ₀,1)-h.Nν_shell)
    χ₀_core = view(χ₀,ind_core)
    λ[:] = -s*sum(χ,dims=[2])[:,1] ./ χ₀_core .+ s
    χ_core = sum(χ)/β^2
    λ[:] = (λ .- U*χ₀_asym)/(1-s*U*χ₀_asym)
    λ_s = -sum((U .* λ .- s*U) .* χ₀_core)/β^2
    χ_out = (χ_core + χ₀_asym*(1+2*λ_s+s*U*χ₀_asym))/(1-U^2 * χ₀_asym^2)
    return χ_out
end


"""
    calc_λ0_impr(type::Symbol, ωgrid::AbstractVector{Int},
                 F::AbstractArray{ComplexF64,3}, χ₀::AbstractArray{ComplexF64,3}, 
                 χ₀_asym::Array{ComplexF64,2}, γ::AbstractArray{ComplexF64,2}, 
                 χ::AbstractArray{ComplexF64,1},
                 U::Float64, β::Float64, h::BSE_Asym_Helper)
    
Calculates improved version of `λ₀ = χ₀ ⋆ F`.
TODO: finish documenation and tests.
"""
function calc_λ0_impr(type::Symbol, ωgrid::AbstractVector{Int},
                 F::AbstractArray{ComplexF64,3}, χ₀::AbstractArray{ComplexF64,3}, 
                 χ₀_asym::Array{ComplexF64,2}, γ::AbstractArray{ComplexF64,2}, 
                 χ::AbstractArray{ComplexF64,1},
                 U::Float64, β::Float64, h::BSE_Asym_Helper)
    s = (type == :ch) ? -1 : +1
    ind_core = (h.Nν_shell+1):(size(χ₀,2)-h.Nν_shell)
    Nq = size(χ₀,1)
    Nν = length(ind_core)
    Nω = size(χ₀,3)
    λasym = Array{ComplexF64,1}(undef, Nν)
    λcore = Array{ComplexF64,1}(undef, Nν)
    res = Array{ComplexF64,3}(undef, Nq, Nν, Nω)

    for (ωi,ωn) in enumerate(ωgrid)
        λasym = -(view(γ,:,ωi) .* (1 .+ s*U .* χ[ωi]) ) .+ 1
        for qi in 1:Nq
            λcore[:] = [s*dot(view(χ₀,qi,ind_core,ωi), view(F,νi,:,ωi))/(β^2) for νi in 1:size(F,1)]
            F_diag!(type, ωn, U, β, χ₀[qi,:,ωi], h)
            res[qi,:,ωi] = λcore + χ₀_asym[qi,ωi].*U.*(λasym .- 1) .+ view(h.diag_asym_buffer, ind_core)
        end
    end
    return res
end
