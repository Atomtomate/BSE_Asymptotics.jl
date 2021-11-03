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
    return I_core, I_corner, union(I_top_bottom, I_left_right)
end

function aux_indices(ind_lst::Vector{CartesianIndex{2}}, ωi::Int, n_iω::Int, n_iν::Int, shift::Int)
    ind1_list = Vector{Int}(undef, length(ind_lst))    # ν' - ν
    ind2_list = Vector{Int}(undef, length(ind_lst))    # ν + ν' + ω
    for (i,ind) in enumerate(ind_lst)
        ind1_list[i] = ωn_to_ωi(ind[1] - ind[2])
        νn, ωn = OneToIndex_to_Freq(ωi, ind[1], n_iω, n_iν, shift)
        νpn, ωn = OneToIndex_to_Freq(ωi, ind[2], n_iω, n_iν, shift)
        ind2_list[i] = ωn_to_ωi(νn + νpn + 1 + ωn)
    end
    return ind1_list, ind2_list
end

function improve_χ!(χsp, χch, χ₀, χsp_asympt, χch_asympt, χpp_asympt, U::Float64, β::Float64, Nν_shell::Int, shift::Int)
    Nν_full = size(χch, 1)
    #TODO: this assumes -nν:nν-1
    Nν_c = size(χch, 1) - 2*Nν_shell
    n_iω   = trunc(Int,size(χch,3)/2)
    n_iν   = trunc(Int,Nν_full/2)

    # internal arrays
    I_core, I_corner, I_ν = shell_indices(Nν_full, Nν_shell)
    I_all = union(I_core, I_corner, I_ν)
    I_aympt = union(I_corner, I_ν)

    Fsp = Array{ComplexF64,2}(undef, Nν_full, Nν_full)
    Fch = similar(Fsp)
    λsp = Array{ComplexF64, 1}(undef, Nν_full)
    λch = similar(λsp)

    for ωi in axes(χch, 3)
        # setup 
        ω_off = -shift*trunc(Int,(ωi-n_iω-1)/2)
        ind1_list_corner, ind2_list_corner = aux_indices(I_corner, ωi, n_iω, n_iν, shift)
        ind1_list_ν, ind2_list_ν = aux_indices(I_ν, ωi, n_iω, n_iν, shift)

        #println(size(Fsp))
        #println(size(χsp))
        #println(size(χ₀))
        for i in I_all
            δ_ννp = Float64(i[1] == i[2])
            #println(i)
            Fsp[i] = - β^2 * (χsp[i,ωi] - δ_ννp*χ₀[i[1],ωi])/(χ₀[i[1],ωi]*χ₀[i[2],ωi])
            Fch[i] = - β^2 * (χch[i,ωi] - δ_ννp*χ₀[i[1],ωi])/(χ₀[i[1],ωi]*χ₀[i[2],ωi])
        end

        # SC
        converged = false
        for i in 1:2
            #TODO: while !converged
            χsp_n = update_χ!(λsp, view(χsp,:,:,ωi), Fsp, view(χ₀,:,ωi), β, I_aympt)
            χch_n = update_χ!(λch, view(χch,:,:,ωi), Fch, view(χ₀,:,ωi), β, I_aympt)
            update_Fsp!(Fsp, λsp, χsp_n, χch_asympt, χsp_asympt, χpp_asympt, U, I_corner, I_ν, ind1_list_corner, ind2_list_corner, ind1_list_ν, ind2_list_ν)
            update_Fch!(Fch, λch, χch_n, χch_asympt, χsp_asympt, χpp_asympt, U, I_corner, I_ν, ind1_list_corner, ind2_list_corner, ind1_list_ν, ind2_list_ν)
            #TODO: check convergence
            converged = true
        end
    end
    return 
end

function update_Fsp!(F, λ, χ, χch_asympt, χsp_asympt, χpp_asympt, U::Float64, I_corner, I_ν, ind1_list_corner, ind2_list_corner, ind1_list_ν, ind2_list_ν)
    for i in 1:length(I_corner)
        i1 = I_corner[i]
        i2 = ind1_list_corner[i]
        i3 = ind2_list_corner[i]
        F[i1] = U + (U^2/2)*χch_asympt[i2] + 3*(U^2/2)*χsp_asympt[i2] - (U^2)*χpp_asympt[i3] + U*λ[i1[1]] + U*λ[i1[2]] + (U^2)*χ    
    end
    for i in 1:length(I_ν)
        i1 = I_ν[i]
        i2 = ind1_list_ν[i]
        i3 = ind2_list_ν[i]
        F[i1] = U + (U^2/2)*χch_asympt[i2] + 3*(U^2/2)*χsp_asympt[i2] - (U^2)*χpp_asympt[i3] + U*λ[i1[1]] + U*λ[i1[2]] 
    end
end

function update_Fch!(F, λ, χ, χch_asympt, χsp_asympt, χpp_asympt, U::Float64, I_corner, I_ν, ind1_list_corner, ind2_list_corner, ind1_list_ν, ind2_list_ν)
    for i in 1:length(I_corner)
        i1 = I_corner[i]
        i2 = ind1_list_corner[i]
        i3 = ind2_list_corner[i]
        F[i1] = -U + (U^2/2)*χch_asympt[i2] - (U^2/2)*χsp_asympt[i2] + (U^2)*χpp_asympt[i3] + U*λ[i1[1]] + U*λ[i1[2]] + (U^2)*χ    
    end
    for i in 1:length(I_ν)
        i1 = I_ν[i]
        i2 = ind1_list_ν[i]
        i3 = ind2_list_ν[i]
        F[i1] = -U + (U^2/2)*χch_asympt[i2] - (U^2/2)*χsp_asympt[i2] + (U^2)*χpp_asympt[i3] + U*λ[i1[1]] + U*λ[i1[2]] 
    end
end

function update_χ!(λ, χ::AbstractArray{ComplexF64,2}, F::AbstractArray{ComplexF64,2}, χ₀::AbstractArray{ComplexF64,1}, β::Float64, indices)
    for i in indices
        (i[1] == i[2]) && (χ[i[1],i[1]] = χ₀[i[1]])
        χ[i] -= χ₀[i[1]]*F[i]*χ₀[i[2]]/(β^2)
    end
    #TODO: absorb into sum
    for νk in axes(F, 3)
        λ[νk]  = sum(χ[:,νk]) / (χ₀[νk]) - 1
        #!(λ[νk] ≈ sum(χ[νk,:]) / (χ₀[νk]) - 1) && error("ν sums not equal: $(λ[νk]) vs $(sum(χ[νk,:]))")
    end
    χs = sum(χ)/β
    return χs
end
