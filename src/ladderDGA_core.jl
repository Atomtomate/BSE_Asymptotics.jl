#TODO: this has to be moved to LadderDGA_core


@inline function OneToIndex_to_Freq(ωi::Int, νi::Int, n_iω::Int, n_iν::Int, shift::Int)
    ωn = ωi-n_iω-1
    νn = (νi-n_iν-1) - shift*trunc(Int,ωn/2)
    return ωn, νn
end

function calc_bubble(Gνω, kG::ReducedKGrid, β::Float64, n_iω::Int, n_iν::Int, shift::Int)
    χ₀ = Array{ComplexF64,3}(undef, 1, 2*n_iν, 2*n_iω+1)
    for ωi in axes(χ₀,3), νi in axes(χ₀,2)
        ωn, νn = OneToIndex_to_Freq(ωi, νi, n_iω, n_iν, shift)
        χ₀[1,νi,ωi] = Gνω[νn]*Gνω[νn+ωn]
        χ₀[1,νi,ωi] .*= -β
    end
    return bubble
end
