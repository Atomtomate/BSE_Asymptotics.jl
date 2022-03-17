"""
Computes asymptotically improved versions of the physical susceptibility χ^{ω}_{r}
and the γ^{ν,ω}_{r} vertex from χ^{ν,νp,ω}_{r}.

Input: χ^{ν,νp,ω}_{r}
Output: χ^{ω}_{r}, λ^{ν,ω}_{r} oder γ^{ν,ω}_{r}

"""
module BSE_SC

export setup, read_gImp, BSE_SC_Helper, BSE_Asym_Helper
export χ₀_shell_sum_core, χ₀_shell_sum, calc_χ₀, improve_χ!, calc_χλ_impr, calc_λ0_impr

using LinearAlgebra
using JLD2
using OffsetArrays
using TimerOutputs

include("ladderDGA_core.jl")
include("helpers.jl")
include("IO.jl")
include("build_chi_asympt.jl")
include("dbg_tools.jl")

function __init__()
    global to = TimerOutput()
end

end
