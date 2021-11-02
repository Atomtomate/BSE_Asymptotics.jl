"""
    =========-===========
    ∥     N_shell       |
    ∥   ___________     |
    ∥   |          |    |
    ∥   |..N_core..|    |
    ∥   |          |    |
    ∥   |__________|    |
    ∥                   |
    =====================

"""
function shell_indices(Nν_full::Int, Nν_shell::Int)
    ind_inner = (Nν_shell+1):(Nν_full-Nν_shell)
    ind_outer = union(1:(Nν_shell),(Nν_full-Nν_shell+1):Nν_full)
    I_top_bottom = [CartesianIndex(νi,νpi) for νi in ind_outer for νpi in ind_inner]
    I_left_right = [CartesianIndex(νi,νpi) for νi in ind_inner for νpi in ind_outer]
    I_corner = [CartesianIndex(νi,νpi) for νi in ind_outer for νpi in ind_outer]
    I_core = [CartesianIndex(νi,νpi) for νi in ind_inner for νpi in ind_inner]
    return I_core, I_corner, I_top_bottom, I_left_right
end

function build_χγ(χsp, χch, χ₀, χch_asympt, χsp_asympt, χpp_asympt, U::Float64, β::Float64, Nν_shell::Int, shift::Int)
    Nν_full = size(χ_ch, 1)
    #TODO: this assumes -nν:nν-1
    Nν_c = size(χ_ch, 1) - 2*Nν_shell
    Nω   = trunc(Int,size(χ_ch,3)/2)

    # internal arrays
    I_core, I_corner, I_ν, I_νp = shell_indices(Nν_full, Nν_shell)
    Fsp = Array{ComplexF64,2}(undef, Nν_full, Nν_full)
    Fch = similar(Fsp)
    χwork = similar(Fch)
    λsp = Array{ComplexF64, 1}(undef, Nν_full)
    λch = similar(λsp)

    for ωi in axes(χ_ch, 3)
        ω_off = shift*trunc(Int,(ωi-Nω-1)/2)
        converged = false
        for i in I_core
            δ_ννp = Float64(i[1] == i[2])
            Fsp[i] = - β^2 * (χsp[i,ωi] - δ_ννp*χ₀[i[1],ωi])/(χ₀[i[1],ωi]*χ₀[i[2],ωi])
            Fch[i] = - β^2 * (χch[i,ωi] - δ_ννp*χ₀[i[1],ωi])/(χ₀[i[1],ωi]*χ₀[i[2],ωi])
        end
        for i in 1:2
        #TODO: while !converged
            χsp_n = calc_χλ!(λsp, χwork, Fsp, view(χ₀,1,:,ωi), β)
            χch_n = calc_χλ!(λch, χwork, Fch, view(χ₀,1,:,ωi), β)
            update_Fsp!(Fsp, λsp,χch_asympt, χsp_asympt, χpp_asympt, U, ωi, Nω, shift, I_corners, I_ν, I_νp)
            #TODO: check convergence
            converged = true
        end
    end
    return Fsp
end

function update_Fsp!(F, λ, χch_asympt, χsp_asympt, χpp_asympt, U::Float64, ωi::Int, Nω, shift, I_corners, I_ν, I_νp)
    ωn = ωi-Nω-1
    shift_ind = - shift*trunc(Int,ωn/2)
    #TODO: construct I2_side, I3_side (nu-nup, nu+nup+omega)
    for i in I_corners
        ind1 = (i[2]-i[1]) + 1
        ind1 = ind1 < 0 ? -ind1 : ind1

    #ωn = ωi-n_iω-1
    #νn = (νi-n_iν-1) - shift*trunc(Int,ωn/2)
    #νpn = (νpi-n_iν-1) - shift*trunc(Int,ωn/2)
    #νi + shift*trunc(Int,ωn/2) + n_iν + 1 = νi
    #νpn + shift*trunc(Int,ωn/2) + n_iν + 1 = νpi
        ind2 = (i[1]+i[2]) + 1
        ind2 = ind1 < 0 ? -ind1 : ind1

        F[i] = U + (U^2/2)*χch_asympt[i[2]-i[1]] + 3*(U^2/2)*χsp_asympt[i[2]-i[1]] - (U^2)*χpp_asympt[i[1]+i[2]+ωi+2*shift_ind]
    end
    for i in I_ν
        F[i] = U + (U^2/2)*χch_asympt[i[2]-i[1]] + 3*(U^2/2)*χsp_asympt[i[2]-i[1]] - (U^2)*χpp_asympt[i[1]+i[2]+ωi+2*shift_ind]
    end
    for i in I_νp
        F[i] = U + (U^2/2)*χch_asympt[i[2]-i[1]] + 3*(U^2/2)*χsp_asympt[i[2]-i[1]] - (U^2)*χpp_asympt[i[1]+i[2]+ωi+2*shift_ind]
    end
end

function calc_χλ!(λ, χwork::AbstractArray{ComplexF64,2}, F, χ₀, β::Float64)
    fill!(χwork, zero(eltype(χwork)))
    for νi in axes(χwork,1)
        χwork[νi,νi] = χ₀[νi]
        for νpi in axes(χwork,2)
            χwork[νi,νpi] -= χ₀[νi]*F[νi,νpi]*χ₀[νpi]/(β^2)
        end
    end
    #TODO: absorb into sum
    for νk in axes(F, 3)
        λνω[νk]  = sum(χwork[:,νk]) / (χ₀[νk]) - 1
        !(λνω[νk] ≈ sum(χwork[νk,:]) / (χ₀[νk]) - 1) && error("ν sums not equal")
    end
    χ = sum(χwork)/β
    return χ
end
