#TODO: index struct
"""
    improve_χ!(type::Symbol, ωi::Int, χr, χ₀, U, β, h::BSE_SC_Helper; 
               Nit=200, atol=1e-9)

Improves asymptotics of `χr`, given a channel `type = :m` or `:d`, using Eq. 12a 12b 
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
    f = if type == :m
        update_Fsp!
    elseif type == :d
        update_Fch!
    else
        error("Unkown channel. Only m/d implemented")
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
        h.Fr[i1] = -U + (U^2)*χ + U*h.λr[i1[1]]  + U*h.λr[i1[2]] + (U^2/2)*h.χch_asympt[i2] - (U^2/2)*h.χsp_asympt[i2] + (U^2)*h.χpp_asympt[i3] 
    end
end

function update_Fch!(χ::ComplexF64, U::Float64, ωn::Int, h::BSE_SC_Helper)
    i1_l = h.ind1_list
    i2_l = view(h.ind2_list, :, ωn)
    for i in 1:length(i1_l)
        i1 = h.I_asympt[i]
        i2 = i1_l[i]
        i3 = i2_l[i]
        h.Fr[i1] = U + (U^2)*χ  - U*h.λr[i1[1]] - U*h.λr[i1[2]] + (U^2/2)*h.χch_asympt[i2] + 3*(U^2/2)*h.χsp_asympt[i2] - (U^2)*h.χpp_asympt[i3]
    end
end

function update_χ!(λ::AbstractArray{ComplexF64,1}, χ::AbstractArray{ComplexF64,2}, 
        F::AbstractArray{ComplexF64,2}, χ₀::AbstractArray{ComplexF64,1}, 
        β::Float64, indices::AbstractArray{CartesianIndex{2},1})
    for i in indices
        χ[i] = i[1] == i[2] ? χ₀[i[1]] : 0
        χ[i] -= χ₀[i[1]]*F[i]*χ₀[i[2]]/(β^2)
    end
    for νk in 1:size(χ,2)
        λ[νk] = (sum(view(χ,:,νk)) / (-χ₀[νk])) + 1
    end
    return sum(χ)/(β^2)
end

function F_diag!(diag_asym_buffer::Vector{ComplexF64}, qi::Int, ωi::Int,  ωn::Int, χ₀::Array{ComplexF64,3},
                        buffer::OffsetMatrix{ComplexF64}, h::BSE_Asym_Helper)

    for i in eachindex(h.block_i)
        ii = h.block_i[i]
        diag_asym_buffer[ii] = 0
        for j in h.block_slices[i]
            diag_asym_buffer[ii] += buffer[j, ωn]*(-χ₀[qi,h.ind1_list[j],ωi])
        end
    end
end

function F_diag!(diag_asym_buffer::Vector{ComplexF64}, qi::Int, ωi::Int,  ωn::Int, χ₀::Array{ComplexF64,3},
                        buffer::Matrix{ComplexF64}, h::BSE_Asym_Helper)

    for i in eachindex(h.block_i)
        ii = h.block_i[i]
        diag_asym_buffer[ii] = 0
        for j in h.block_slices[i]
            diag_asym_buffer[ii] += buffer[j, ωi]*(-χ₀[qi,h.ind1_list[j],ωi])
        end
    end
end

F_diag!(type, ωn, U, β, χ₀, h::BSE_Asym_Helper_Approx2) = nothing


"""
    calc_χλ(type::Symbol, ωn::Int, χ::AbstractArray{ComplexF64,2}, χ₀::AbstractArray{ComplexF64,1}, U::Float64, β::Float64, χ₀_asym::Float64, h::BSE_Asym_Helper)

Calculates the physical susceptibility `χ` and triangular vertex `λ` in a given channel `type=:m` or `type=:d` using knowledge about the asymptotics of the full vertex and tails of the Green's function.
`χ₀_asym` is the `ω` dependent asymptotic tail of `χ₀` and can be calculated with  [`χ₀_shell_sum`](@ref).

TODO: refactor code duplications
"""
function calc_χλ_impr(type::Symbol, qi::Int, ωi::Int, ωn::Int, χ::Array{ComplexF64,2}, χ₀::Array{ComplexF64,3}, 
                 U::Float64, β::Float64, χ₀_asym::ComplexF64, h::HT) where  HT <: BSE_Asym_Helpers
    λ = Array{eltype(χ),1}(undef, size(χ, 1))
    χ_out = calc_χλ_impr!(λ, h.diag_asym_buffer, type, qi, ωi, ωn, χ, χ₀, U, β, χ₀_asym, h)
    return χ_out, λ
end

function calc_χλ_impr!(λ::Vector{ComplexF64}, diag_asym_buffer::Vector{ComplexF64}, type::Symbol, qi::Int, ωi::Int, ωn::Int, χ::Array{ComplexF64,2}, χ₀::Array{ComplexF64,3}, 
                 U::Float64, β::Float64, χ₀_asym::ComplexF64, h::BSE_Asym_Helper)::ComplexF64
    s = type == :d ? -1 : 1
    F_diag!(diag_asym_buffer, qi, ωi, ωn, χ₀, type == :d ? h.buffer_d : h.buffer_m, h)

    core_offset::Int = h.Nν_shell
    norm::ComplexF64 = 1 / (1-s*U*χ₀_asym)
    λ_s::ComplexF64 = 0.0
    χ_core::ComplexF64 = 0.0
    for νi in axes(λ, 1)
        λ_core::ComplexF64 = sum(view(χ,νi, :))
        χ_core += λ_core
        χ₀_val::ComplexF64 = χ₀[qi,core_offset+νi,ωi]
        @inbounds λ[νi] = (-s*λ_core / χ₀_val + s - s * diag_asym_buffer[core_offset+νi] - U*χ₀_asym) * norm
        @inbounds λ_s -= (U * λ[νi] - s*U) * χ₀_val
    end
    λ_s /= β^2
    χ_core /= β^2
    diag_asym_s::ComplexF64 = -dot(diag_asym_buffer, view(χ₀,qi,:,ωi))/β^2
    χ_out::ComplexF64 = (χ_core + χ₀_asym*(1+2*λ_s+s*U*χ₀_asym) - diag_asym_s)/(1-U^2 * χ₀_asym^2)

    return χ_out
end


"""
    calc_λ0_impr(type::Symbol, ωgrid::AbstractVector{Int},
                 F::AbstractArray{ComplexF64,3}, χ₀::AbstractArray{ComplexF64,3}, 
                 χ₀_asym::Array{ComplexF64,2}, γ::AbstractArray{ComplexF64,2}, 
                 χ::AbstractArray{ComplexF64/Float64,1},
                 U::Float64, β::Float64, h::BSE_Asym_Helper; 
                 diag_zero::Bool=true, use_threads::Bool=true)
    
Calculates improved version of `λ₀ = χ₀ ⋆ F`.
TODO: finish documentation and tests.
"""
function calc_λ0_impr(type::Symbol, ωgrid::AbstractVector{Int},
                 F::AbstractArray{ComplexF64,3}, χ₀::AbstractArray{ComplexF64,3}, 
                 χ₀_asym::Array{ComplexF64,2}, γ::AbstractArray{ComplexF64,2}, 
                 χ::AbstractArray{T,1},
                 U::Float64, β::Float64, h; diag_zero::Bool=true, use_threads::Bool=true) where T <: Union{ComplexF64,Float64}
    @assert type in [:m, :d] "type must be :m or :d"

    s = (type == :d) ? -1 : +1
    ind_core = (h.Nν_shell+1):(size(χ₀,2)-h.Nν_shell)
    Nq = size(χ₀,1)
    Nν = length(ind_core)
    Nω = size(χ₀,3)
    core_offset::Int = h.Nν_shell
    

    res = Array{ComplexF64,3}(undef, Nq, Nν, Nω)

    diag_term = if !diag_zero && hasfield(typeof(h), :diag_asym_buffer)
        println(stderr, "DBG: using diagonal terms in λ₀")
    else
        println(stderr, "DBG: NOT using diagonal terms in λ₀")
    end

    
    if use_threads
        NT = Threads.nthreads()
        diag_buffers = [zeros(ComplexF64,size(h.diag_asym_buffer)) for ti in 1:NT]
        Threads.@threads for ωi in eachindex(ωgrid)
            ωn = ωgrid[ωi]
            for νi in axes(F,1)
                for qi in 1:Nq
                    if !diag_zero && νi == 1 && hasfield(typeof(h), :diag_asym_buffer)
                        F_diag!(diag_buffers[Threads.threadid()], qi, ωi, ωn, χ₀, type == :d ? h.buffer_d : h.buffer_m, h)
                    end
                    λcore = s*dot(view(χ₀,qi,ind_core,ωi), view(F,νi,:,ωi))/(β^2) 
                    λasym = -(γ[νi,ωi] * (1 + s*U * χ[ωi]) ) + 1
                    res[qi,νi,ωi] = λcore + χ₀_asym[qi,ωi] * U *(λasym - 1) + diag_buffers[Threads.threadid()][core_offset+νi]
                end
            end
        end
    else
        NT = Threads.nthreads()
        diag_buffer = zeros(ComplexF64,size(h.diag_asym_buffer))
        for ωi in eachindex(ωgrid)
            ωn = ωgrid[ωi]
            for νi in axes(F,1)
                for qi in 1:Nq
                    if !diag_zero && νi == 1 && hasfield(typeof(h), :diag_asym_buffer)
                        F_diag!(diag_buffer, qi, ωi, ωn, χ₀, type == :d ? h.buffer_d : h.buffer_m, h)
                    end
                    λcore = s*dot(view(χ₀,qi,ind_core,ωi), view(F,νi,:,ωi))/(β^2) 
                    λasym = -(γ[νi,ωi] * (1 + s*U * χ[ωi]) ) + 1
                    res[qi,νi,ωi] = λcore + χ₀_asym[qi,ωi] * U *(λasym - 1) + diag_buffer[core_offset+νi]
                end
            end
        end
    end
    return res
end
