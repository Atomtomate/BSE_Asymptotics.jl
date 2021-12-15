function gen_synth(Nf::Int, Nb::Int, Nbor::Int, shift::Int)
    gi = [1,2,3]
    χ₀ = Float64.(collect(reshape(100*(1:((2*(Nf+Nbor))*(2*Nb+1))),2*(Nf+Nbor),2*Nb+1)))
    χsp = Float64.(collect(reshape(1:((2*(Nf+Nbor))^2*(2*Nb+1)),2*(Nf+Nbor),2*(Nf+Nbor),2*Nb+1)))
    χch = Float64.(collect(reshape(1:((2*(Nf+Nbor))^2*(2*Nb+1)),2*(Nf+Nbor),2*(Nf+Nbor),2*Nb+1)))
    χsp_a = Float64.(collect(1000:(1000+Nf+100)))
    χch_a = Float64.(collect(1000:(1000+Nf+100)))
    χpp_a = Float64.(collect(1000:(1000+Nf+100)))
    U = 1.0
    β = 1.0
    gi, χ₀, χsp, χch, χsp_a, χch_a, χpp_a, U, β, shift 
end

function update_Fsp!(χ::ComplexF64, U::Float64, ωi::Int, h)
    i1_l = h.ind1_list
    i2_l = view(h.ind2_list, :, ωi)
    for i in 1:length(h.I_asympt)
        i1 = h.I_asympt[i]
        i2 = i1_l[i]
        i3 = i2_l[i]
        h.Fr[i1] = -U + U*h.λr[i1[1]]  + U*h.λr[i1[2]] + (U^2)*χ 
        #+ (U^2/2)*h.χch_asympt[i2] - (U^2/2)*h.χsp_asympt[i2] + (U^2)*h.χpp_asympt[i3]  + 
    end
end

function update_Fch!(χ::ComplexF64, U::Float64, ωi::Int, h)
    i1_l = h.ind1_list
    i2_l = view(h.ind2_list, :, ωi)
    for i in 1:length(h.I_asympt)
        i1 = h.I_asympt[i]
        i2 = i1_l[i]
        i3 = i2_l[i]
        h.Fr[i1] = U + (U^2)*χ  - U*h.λr[i1[1]] - U*h.λr[i1[2]] #+
              #(U^2/2)*h.χch_asympt[i2] + 3*(U^2/2)*h.χsp_asympt[i2] - (U^2)*h.χpp_asympt[i3]
    end
end

function improve_χ_trace!(type::Symbol, ωi::Int, χr::AbstractArray{ComplexF64,2}, χ₀::AbstractArray{ComplexF64,1}, 
                U::Float64, β::Float64, h; Nit=200, atol=1e-9)
    f = if type == :sp
        update_Fsp!
    elseif type == :ch
        update_Fch!
    else
        error("Unkown channel. Only sp/ch implemented")
    end

    Fr_trace = []
    χr_trace = []
    χlocr_trace = []
    λr_trace = []
    fill!(h.Fr, 0.0)
    for i in h.I_core
        δ_ννp = Float64(i[1] == i[2])
        h.Fr[i] = - β^2 * (χr[i] - δ_ννp*χ₀[i[1]])/(χ₀[i[1]]*χ₀[i[2]])
    end
    χr_old = 0.0
    converged = false
    i = 0
    push!(Fr_trace, deepcopy(h.Fr))
    push!(χr_trace, deepcopy(χr))
    push!(χlocr_trace, deepcopy(χr_old))
    push!(λr_trace, deepcopy(h.λr))
    while !converged && (i < Nit)
        χr_n = update_χ!(h.λr, χr, h.Fr, χ₀, β, h.I_asympt)
        f(χr_n, U, ωi, h)
        if (abs(χr_old - χr_n) < atol)
            converged = true
        else
            i += 1
            χr_old = χr_n
        end
        if i == 1 || i == 3 || (i == Nit || converged)
            push!(Fr_trace, deepcopy(h.Fr))
            push!(χr_trace, deepcopy(χr))
            push!(χlocr_trace, deepcopy(χr_old))
            push!(λr_trace, deepcopy(h.λr))
    end
    end
    return i, Fr_trace, χr_trace, χlocr_trace, λr_trace
end


function improve_χλ_direct(ωn::Int, χsp::AbstractArray{ComplexF64,2}, χch::AbstractArray{ComplexF64,2}, χ₀::AbstractArray{ComplexF64,1}, U::Float64, β::Float64, bs, h::BSE_Asym_Helper)
    i1_l = h.ind1_list
    i2_l = view(h.ind2_list, :, ωn)
    ind_core = (h.Nν_shell+1):(length(χ₀)-h.Nν_shell)
    χ₀_core = view(χ₀,ind_core)
    λsp_core = -sum(χsp,dims=[2])[:,1] ./ χ₀_core .+ 1
    λch_core =  sum(χch,dims=[2])[:,1] ./ χ₀_core .- 1
    χsp_core = sum(χsp) /β^2
    χch_core = sum(χch) /β^2
    diag_asym_sp = zeros(eltype(h.χch_asympt), length(χ₀))
    diag_asym_ch = zeros(eltype(h.χch_asympt), length(χ₀))
    #TODO this is inefficient. only loop over necessary indices, only one direction needs to be extended
    for i in 1:length(i1_l)
        i1 = h.I_asympt[i]
        i2 = i1_l[i]
        i3 = i2_l[i]
        diag_asym_sp[i1[1]] += ((U^2/2)*h.χch_asympt[i2] -   (U^2/2)*h.χsp_asympt[i2] + (U^2)*h.χpp_asympt[i3])*(-χ₀[i1[2]])/β^2
        diag_asym_ch[i1[1]] += ((U^2/2)*h.χch_asympt[i2] + 3*(U^2/2)*h.χsp_asympt[i2] - (U^2)*h.χpp_asympt[i3])*(-χ₀[i1[2]])/β^2
    end
    λsp = (λsp_core .- view(diag_asym_sp,ind_core) .+ U*bs)/(1+U*bs)
    λch = (λch_core .+ view(diag_asym_ch,ind_core) .+ U*bs)/(1-U*bs)
    λsp_s = -sum((U .* λsp .- U) .* χ₀_core)/β^2
    λch_s = -sum((U .* λch .+ U) .* χ₀_core)/β^2
    diag_asym_sp_s = -sum(diag_asym_sp .* χ₀)/β^2
    diag_asym_ch_s = -sum(diag_asym_ch .* χ₀)/β^2
    χsp = (χsp_core - bs*(1+2*λsp_s-U*bs) - diag_asym_sp_s)/(1-U^2 * bs^2)
    χch = (χch_core - bs*(1+2*λch_s+U*bs) - diag_asym_ch_s)/(1-U^2 * bs^2)
    return χsp, χch, λsp, λch, λsp_core, λch_core
end

function calc_χ₀_shell_sum_old(β::Float64, c2::Float64, c3::Float64, n_iω::Int, n_iν::Int, shift::Int)
    χ₀_shell_sum = zeros(ComplexF64, 2*n_iω+1)
    iω_grid = iω_array(β, -n_iω:n_iω)
    for (ωi,ωn) in enumerate(-n_iω:n_iω)
        si = shift*trunc(Int,ωn/2)
        χ₀_shell_sum[ωi] = 0
        iν_ind_grid = (-n_iν-n_iω-si):(n_iν+n_iω-si-1)
        iν_grid = OffsetArray(iν_array(β, iν_ind_grid), iν_ind_grid)
        for νn in (-n_iν-si):(n_iν-si-1)
            χ₀_shell_sum[ωi] += ( -1/(iν_grid[νn]^1 * iν_grid[νn+ωn]^1)
                                 +c2/(iν_grid[νn]^2 * iν_grid[νn+ωn]^1)
                                 +c2/(iν_grid[νn]^1 * iν_grid[νn+ωn]^2)
                                 -c3/(iν_grid[νn]^1 * iν_grid[νn+ωn]^3)
                                 -c3/(iν_grid[νn]^3 * iν_grid[νn+ωn]^1)
                                 -c2^2/(iν_grid[νn]^2 * iν_grid[νn+ωn]^2))
        end
        χ₀_shell_sum[ωi] /= β
        χ₀_shell_sum[ωi] += ωn == 0  ? -(β/4-c3*β^3/24-c2^2*β^3/48) : ((c2^2-c3)*β/(2*(2*ωn*π/β)^2))
    end
    return χ₀_shell_sum
end
