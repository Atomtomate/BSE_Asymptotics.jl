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

function improve_χ_trace!(type::Symbol, ωi::Int, χr::AbstractArray{ComplexF64,2}, χ₀::AbstractArray{ComplexF64,1}, 
                U::Float64, β::Float64, shift::Int, h::BSE_SC_Helper; Nit=200, atol=1e-9)
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
        #println("$i: $(real(χr_old)) -> $(real(χr_n))")
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
