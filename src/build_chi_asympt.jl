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
    I_top_bottom = sort([CartesianIndex(νi,νpi) for νi in ind_outer for νpi in ind_inner])
    I_left_right = sort([CartesianIndex(νi,νpi) for νi in ind_inner for νpi in ind_outer])
    I_corner = sort([CartesianIndex(νi,νpi) for νi in ind_outer for νpi in ind_outer])
    I_core = sort([CartesianIndex(νi,νpi) for νi in ind_inner for νpi in ind_inner])
    return I_core, I_corner, I_top_bottom, I_left_right
end

function χ₀_indices(ind_lst::Vector{CartesianIndex{2}}, ωi::Int, n_iω::Int, n_iν::Int, shift::Int)
    ind1_list = Vector{Int}(undef, length(ind_lst))    # ν
    ind2_list = Vector{Int}(undef, length(ind_lst))    # ν'
    for (i,ind) in enumerate(ind_lst)
        ωn, νn, = OneToIndex_to_Freq(ωi, ind[1], n_iω, n_iν, shift)
        ωn, νpn = OneToIndex_to_Freq(ωi, ind[2], n_iω, n_iν, shift)
        ind1_list[i] = ωn_to_ωi(νn)
        ind2_list[i] = ωn_to_ωi(νpn)
    end
    return ind1_list, ind2_list
end

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

function improve_χ!(χsp, χch, χ₀, χsp_asympt, χch_asympt, χpp_asympt, U::Float64, β::Float64, Nν_shell::Int, shift::Int; active_terms=[1,1,1,1])
    Nν_full = size(χch, 1)
    #TODO: this assumes -nν:nν-1
    Nν_c = size(χch, 1) - 2*Nν_shell
    n_iω   = trunc(Int,size(χch,3)/2)
    n_iν   = trunc(Int,Nν_full/2)

    # internal arrays
    I_core, I_corner, I_t, I_r = shell_indices(Nν_full, Nν_shell)
    I_all = sort(union(I_core, I_corner, I_r, I_t))
    I_aympt = sort(union(I_corner, I_r, I_t))

    Fsp_trace = []
    Fch_trace = []
    Fsp = Array{ComplexF64,2}(undef, Nν_full, Nν_full)
    Fch = similar(Fsp)
    λsp = Array{ComplexF64, 1}(undef, Nν_full)
    λch = similar(λsp)

    for ωi in axes(χch,3) #(n_iω+1):(n_iω+1)#
        print(ωi )
        flush(stdout)
        # setup 
        #TODO: test w_offset here
        ω_off = shift*trunc(Int,(ωi-n_iω-1)/2)
        ind1_list_corner, ind2_list_corner = aux_indices(I_corner, ωi, n_iω, n_iν, shift)
        ind1_list_r, ind2_list_r = aux_indices(I_r, ωi, n_iω, n_iν, shift)
        ind1_list_t, ind2_list_t = aux_indices(I_t, ωi, n_iω, n_iν, shift)

        #println(size(χ₀))
        fill!(Fsp, 0.0)
        fill!(Fch, 0.0)
        for i in I_core
            δ_ννp = Float64(i[1] == i[2])
            Fsp[i] = - β^2 * (χsp[i,ωi] - δ_ννp*χ₀[i[1],ωi])/(χ₀[i[1],ωi]*χ₀[i[2],ωi])
            Fch[i] = - β^2 * (χch[i,ωi] - δ_ννp*χ₀[i[1],ωi])/(χ₀[i[1],ωi]*χ₀[i[2],ωi])
        end
        # SC
        Nit = 200
        χsp_old = 0.0
        χch_old = 0.0
        converged = false
        i = 0
        while !converged && (i < Nit)
            (ωi == (n_iω + 1)) && (i in [1,2,3,4,Nit]) && push!(Fsp_trace, deepcopy(Fsp))
            (ωi == (n_iω + 1)) && (i in [1,2,3,4,Nit]) && push!(Fch_trace, deepcopy(Fch))
            #TODO: while !converged
            χsp_n = update_χ!(λsp, view(χsp,:,:,ωi), Fsp, view(χ₀,:,ωi), β, I_aympt)
            χch_n = update_χ!(λch, view(χch,:,:,ωi), Fch, view(χ₀,:,ωi), β, I_aympt)
            update_Fsp!(Fsp, λsp, χsp_n, χch_asympt, χsp_asympt, χpp_asympt, U, I_corner, I_r, I_t, ind1_list_corner, ind2_list_corner)
            update_Fch!(Fch, λch, χch_n, χch_asympt, χsp_asympt, χpp_asympt, U, I_corner, I_r, I_t, ind1_list_corner, ind2_list_corner) 
            ((ωi == n_iω+1) && i < 50) && println("$(real(χsp_n)), $(minimum(real.(λsp)))<->$(maximum(real.(λsp))) : $(real(χch_n)),$(minimum(real.(λch)))<->$(maximum(real.(λch)))")
            #TODO: check convergence
            if (abs(χch_old - χch_n) < atol) && (abs(χsp_old - χsp_n) < atol)
                converged = true
            else
                i += 1
                χch_old = χch_n
                χsp_old = χsp_n
            end
        end
    end
    return Fsp_trace, Fch_trace
end

function update_Fsp!(F, λ, χ, χch_asympt, χsp_asympt, χpp_asympt, U::Float64, I_corner, I_r, I_t, ind1_list_corner, ind2_list_corner)
    for i in 1:length(I_corner)
        i1 = I_corner[i]
        i2 = ind1_list_corner[i]
        i3 = ind2_list_corner[i]
        F[i1] = -U + (U^2/2)*χch_asympt[i2] - (U^2/2)*χsp_asympt[i2] + (U^2)*χpp_asympt[i3] - (U^2)*χ # + U*λ[i1[1]]  + U*λ[i1[2]]
    end
    for i in 1:length(I_r)
        i1 = I_r[i]
        F[i1] = -U  + U*λ[i1[1]]
    end
    for i in 1:length(I_t)
        i1 = I_t[i]
        F[i1] = -U  + U*λ[i1[2]]
    end
end

function update_Fch!(F, λ, χ, χch_asympt, χsp_asympt, χpp_asympt, U::Float64, I_corner, I_r, I_t, ind1_list_corner, ind2_list_corner)
    for i in 1:length(I_corner)
        i1 = I_corner[i]
        i2 = ind1_list_corner[i]
        i3 = ind2_list_corner[i]
        F[i1] = U + (U^2/2)*χch_asympt[i2] + 3*(U^2/2)*χsp_asympt[i2] - (U^2)*χpp_asympt[i3] - (U^2)*χ #- U*λ[i1[1]] - U*λ[i1[2]]
    end
    for i in 1:length(I_r)
        i1 = I_r[i]
        F[i1] = U  - U*λ[i1[1]]
    end
    for i in 1:length(I_t)
        i1 = I_t[i]
        F[i1] = U  - U*λ[i1[2]]
    end
end

function update_χ!(λ, χ::AbstractArray, F::AbstractArray, χ₀::AbstractArray, β::Float64, indices)
    for i in indices
        if i[1] == i[2] 
            χ[i] = χ₀[i[1]]
        else
            χ[i] = 0.0
        end
        χ[i] -= χ₀[i[1]]*F[i]*χ₀[i[2]]/(β^2)
    end
    #TODO: absorb into sum
    for νk in 1:size(χ,2)
        λ[νk]  = (sum(χ[:,νk]) / (-χ₀[νk])) + 1
        #λt = sum(χ[νk,:]) / (-χ₀[νk]) + 1
        #!isapprox(λ[νk] , λt, atol=0.00001) && error("ν sums not equal: $(λ[νk]) vs $(λt)")
    end
    χs = sum(χ)/(β^2)
    return χs
end
