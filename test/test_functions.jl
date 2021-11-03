function gen_synth(Nf::Int, Nb::Int, Nbor::Int, shift::Int)
    gi = [1,2,3]
    χ₀ = reshape(100*(1:((2*(Nf+Nbor))^2*(2*Nb+1))),2*(Nf+Nbor),2*(Nf+Nbor),2*Nb+1)
    χsp = reshape(1:((2*(Nf+Nbor))^2*(2*Nb+1)),2*(Nf+Nbor),2*(Nf+Nbor),2*Nb+1)
    χch = reshape(1:((2*(Nf+Nbor))^2*(2*Nb+1)),2*(Nf+Nbor),2*(Nf+Nbor),2*Nb+1)
    χsp_a = collect(1000:(1000+Nf+100))
    χch_a = collect(1000:(1000+Nf+100))
    χpp_a = collect(1000:(1000+Nf+100))
    U = 1
    β = 1
    gi, χ₀, χsp, χch, χsp_a, χch_a, χpp_a, U, β, shift 
end
