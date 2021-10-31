module BSE_SC
"""
Computes asymptotically improved versions of the physical susceptibility χ^{ω}_{r}
and the γ^{ν,ω}_{r} vertex from χ^{ν,νp,ω}_{r}.

Input: χ^{ν,νp,ω}_{r}
Output: χ^{ω}_{r}, λ^{ν,ω}_{r} oder γ^{ν,ω}_{r}

"""

include("ladderDGA_core.jl")
include("LapackWrapper.jl")
include("build_chi_asympt.jl")

function F_d(U::Float64, χ_d, χ_m, χ_pp, λ_d, ω::Int, ν::Int, νp::Int)
    U + (U^2)/2 * χ_d(νp-ν) + 3*(U^2)/2 * χ_m(νp - ν) - U^2 * χ_pp(ν+νp+ω) + U * λ_d(ν,ω) + U * λ_d(νp, ω) + U^2 * χ_d(ω)
end

function F_m(U::Float64, χ_d, χ_m, χ_pp, λ_m, ω::Int, ν::Int, νp::Int)
    - U + (U^2)/2 * χ_d(νp-ν) + (U^2)/2 * χ_m(νp - ν) + U^2 * χ_pp(ν+νp+ω) + U * λ_m(ν,ω) + U * λ_d(νp, ω) + U^2 * χ_d(ω)
end # module
