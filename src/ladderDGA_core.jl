#TODO: this has to be moved to LadderDGA_core

"""
    ωn_to_ωi(ωn::Int)

Translates Matsubara index (i.e. n ∈ [-N,N]) to Julia array index (i.e. i ∈ [1,2N+1]),
assuming symmetry around 0!
"""
@inline ωn_to_ωi(ωn::Int) = (ωn < 0 ? (-ωn+1) : (ωn+1)) 

"""
    OneToIndex_to_Freq(ωi::Int, νi::Int, n_iω::Int, n_iν::Int, shift::Int)

Translates Julia array indices for bosonic and fermionic axis (i.e. ωi ∈ [1,2`n_iω`+1], νi ∈ [1,2`n_iν`]) to 
Matsubara index (i.e. ωn ∈ [-`n_iω`,`n_iω`], νn ∈ [-`n_iν` -`shift`*ωn/2, `n_iν` -`shift`*ωn/2]).
`shift` is used to center the non-asymptotic core of the vertex related functions in the center of the array.
"""
@inline function OneToIndex_to_Freq(ωi::Int, νi::Int, n_iω::Int, n_iν::Int, shift::Int)
    ωn = ωi-n_iω-1
    νn = (νi-n_iν-1) - shift*trunc(Int,ωn/2)
    return ωn, νn
end

"""
    calc_χ₀(Gνω::OffsetArray, β::Float64, n_iω::Int, n_iν::Int, shift::Int)

Calculates the bubble term `χ₀` from a given Green's function `Gνω`. This is mainly used for standalone
testing and should generally be calculated independenty in the method using this package.
"""
function calc_χ₀(Gνω::OffsetArray, β::Float64, n_iω::Int, n_iν::Int, shift::Int)
    χ₀ = Array{ComplexF64,2}(undef, 2*n_iν, 2*n_iω+1)
    for (ωi,ωn) in enumerate(-n_iω:n_iω)
        for (νi,νn) in enumerate((-n_iν - shift*trunc(Int,ωn/2)):(n_iν - shift*trunc(Int,ωn/2) - 1))
            @inbounds χ₀[νi,ωi] = -β*Gνω[νn]*Gνω[νn+ωn]
        end
    end
    return χ₀
end

iν_array(β::Real, grid::AbstractArray{Int64,1}) = ComplexF64[1.0im*((2.0 *el + 1)* π/β) for el in grid]
iω_array(β::Real, grid::AbstractArray{Int64,1}) = ComplexF64[1.0im*((2.0 *el)* π/β) for el in grid]
