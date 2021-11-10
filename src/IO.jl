function read_gImp(inf::String)
    jldopen(inf,"r") do f
        OffsetArray(cat(conj(reverse(f["gImp"])),f["gImp"], dims=1), (-length(f["gImp"])):(length(f["gImp"])-1))
    end
end

function setup(input_file::String, N_shell::Int)
    #TODO: change lDGAPostprocessing to include shift, MF_grid
    gImp = read_gImp(input_file)
    χ₀, χDMFTsp, χDMFTch, χ_sp_asympt, χ_ch_asympt, χ_pp_asympt, U, β, shift = jldopen(input_file, "r") do f
        χ₀ = calc_χ₀(gImp, f["β"], f["grid_nBose"], f["grid_nFermi"]+N_shell, f["grid_shift"])
        χ₀, permutedims(f["χDMFTsp"], (2,3,1)), permutedims(f["χDMFTch"], (2,3,1)), 2*f["χ_sp_asympt"] ./ f["β"]^2, 2*f["χ_ch_asympt"] ./ f["β"]^2, f["χ_pp_asympt"] ./ f["β"]^2, f["U"], f["β"], f["grid_shift"]
    end
    Nν_full = size(χDMFTsp,1) + 2*N_shell
    Nω = size(χDMFTsp,3)
    χ_sp_improved = zeros(eltype(χDMFTsp), Nν_full, Nν_full, Nω)
    χ_ch_improved = zeros(eltype(χDMFTch), Nν_full, Nν_full, Nω)
    χ_sp_improved[(N_shell+1):(end-N_shell),(N_shell+1):(end-N_shell),:] = χDMFTsp
    χ_ch_improved[(N_shell+1):(end-N_shell),(N_shell+1):(end-N_shell),:] = χDMFTch
    #χ_sp_improved[(N_shell+1):(end-N_shell),(N_shell+1):(end-N_shell),:] = χDMFTsp
    #χ_ch_improved[(N_shell+1):(end-N_shell),(N_shell+1):(end-N_shell),:] = χDMFTch
    return gImp, χ₀, χ_sp_improved, χ_ch_improved, χ_sp_asympt, χ_ch_asympt, χ_pp_asympt, U, β, shift
end
