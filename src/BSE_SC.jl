"""
Computes asymptotically improved versions of the physical susceptibility χ^{ω}_{r}
and the γ^{ν,ω}_{r} vertex from χ^{ν,νp,ω}_{r}.

Input: χ^{ν,νp,ω}_{r}
Output: χ^{ω}_{r}, λ^{ν,ω}_{r} oder γ^{ν,ω}_{r}

"""
module BSE_SC

export read_gImp
export calc_χ₀

using LinearAlgebra
using JLD2
using OffsetArrays
using PaddedViews

include("LapackWrapper.jl")
include("ladderDGA_core.jl")
include("IO.jl")
include("build_chi_asympt.jl")


end
