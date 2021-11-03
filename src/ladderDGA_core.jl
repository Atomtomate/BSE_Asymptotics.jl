#TODO: this has to be moved to LadderDGA_core

@inline ωn_to_ωi(ωn::Int) = (ωn < 0 ? (-ωn+1) : (ωn+1)) 

@inline function OneToIndex_to_Freq(ωi::Int, νi::Int, n_iω::Int, n_iν::Int, shift::Int)
    ωn = ωi-n_iω-1
    νn = (νi-n_iν-1) - shift*trunc(Int,ωn/2)
    return ωn, νn
end

function calc_χ₀(Gνω::OffsetArray, β::Float64, n_iω::Int, n_iν::Int, shift::Int)
    #χ₀ = OffsetArray(Array{Number,3}(undef, 1, 4*n_iν+1, 2*n_iω+1), 1:1, -2*n_iν:2*n_iν, -n_iω:n_iω)
    χ₀ = Array{ComplexF64,2}(undef, 2*n_iν, 2*n_iω+1)
    for (ωi,ωn) in enumerate(-n_iω:n_iω)
        for (νi,νn) in enumerate((-n_iν - shift*trunc(Int,ωn/2)):(n_iν - shift*trunc(Int,ωn/2) - 1))
            χ₀[νi,ωi] = -β*Gνω[νn]*Gνω[νn+ωn]
        end
    end
    return χ₀
end
