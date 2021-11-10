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



function improve_χ!(χsp, χch, χ₀, χsp_asympt, χch_asympt, χpp_asympt, U::Float64, β::Float64, Nν_shell::Int, shift::Int)
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
