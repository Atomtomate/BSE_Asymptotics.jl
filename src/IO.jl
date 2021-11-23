"""
    read_gImp(input_file::String)

Read impurity Green's function from JLD2 file at path `input_file`.
"""
function read_gImp(input_file::String)
    jldopen(input_file,"r") do f
        OffsetArray(cat(conj(reverse(f["gImp"])),f["gImp"], dims=1), (-length(f["gImp"])):(length(f["gImp"])-1))
    end
end

"""
    setup(input_file::String, N_shell::Int)

Setup calculation for data from file. This is used for testing from precomputed data.
The file is expected to contain the following fields: `β`, `U`, `gImp`, `grid_nBose`, `grid_nFermi`,
`grid_shift`, `χDMFTsp`, `χDMFTch`, `χ_sp_asympt`, `χ_ch_asympt`, `χ_pp_asympt`.
"""
function setup(input_file::String, N_shell::Int)
    #TODO: change lDGAPostprocessing to include shift, MF_grid
    gImp = read_gImp(input_file)
    χ₀, χDMFTsp, χDMFTch, χ_sp_asympt, χ_ch_asympt, χ_pp_asympt, U, β, shift = jldopen(input_file, "r") do f
        χ₀ = calc_χ₀(gImp, f["β"], f["grid_nBose"], f["grid_nFermi"]+N_shell, f["grid_shift"])
        χ₀, permutedims(f["χDMFTsp"], (2,3,1)), permutedims(f["χDMFTch"], (2,3,1)), f["χ_sp_asympt"] ./ f["β"]^2, f["χ_ch_asympt"] ./ f["β"]^2, f["χ_pp_asympt"] ./ f["β"]^2, f["U"], f["β"], f["grid_shift"]
    end
    Nν_full = size(χDMFTsp,1) + 2*N_shell
    Nω = size(χDMFTsp,3)
    χ_sp_improved = zeros(eltype(χDMFTsp), Nν_full, Nν_full, Nω)
    χ_ch_improved = zeros(eltype(χDMFTch), Nν_full, Nν_full, Nω)
    χ_sp_improved[(N_shell+1):(end-N_shell),(N_shell+1):(end-N_shell),:] = χDMFTsp
    χ_ch_improved[(N_shell+1):(end-N_shell),(N_shell+1):(end-N_shell),:] = χDMFTch

    n_iω   = trunc(Int,size(χDMFTch,3)/2)
    n_iν   = trunc(Int,Nν_full/2)
    helper = BSE_SC_Helper(χ_sp_asympt, χ_ch_asympt, χ_pp_asympt, Nν_full, N_shell, n_iω, n_iν, shift)
    return gImp, χ₀, χ_sp_improved, χ_ch_improved, helper, U, β
end
