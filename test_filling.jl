using DelimitedFiles
using SpecialFunctions

gm_wim_test_file = "/home/julian/Hamburg/out_of_hf_test/gm_wim_test"

function filling_pos(G::Vector, U::Float64, μ::Float64, β::Float64, shell::Float64)::Float64
    sG = sum(G)
    2*(real(sG + conj(sG))/β + 0.5 + μ * shell) / (1 + U * shell)
end
function filling(G::Vector, U::Float64, μ::Float64, β::Float64, shell::Float64)::Float64
    2*(real(sum(G))/β + 0.5 + μ * shell) / (1 + U * shell)
end

function filling(G::Vector, U::Float64, μ::Float64, β::Float64; ll=1)
    N = floor(Int, length(G)/2)
    shell = shell_sum(N, β)
    filling(G, U, μ, β, shell, ll=ll)
end

"""
    shell_sum(N::Int, β::Float64)::Float64

Calculate ``\\frac{1}{\\beta} \\sum_{n \\in \\Omega_\\mathrm{shell}} \\frac{1}{(i \\nu_n)^2}``
"""
function shell_sum(N::Int, β::Float64)::Float64
    (polygamma(1, N + 1/2) - polygamma(1, 1/2 - N)) * β / (4*π^2) + β/4
end

function shell_sum_naive(iν_array::Vector{ComplexF64}, β::Float64)::Float64
    real(sum((iν_array) .^ 2))/β + β/4
end

r   = readdlm(gm_wim_test_file)[1:10000,:]
v   = 1 ./ (1im .* cat(-r[end:-1:1,1],r[:,1], dims=1))
v_h   = 1 ./ (1im .* r[:,1])
gt  = cat(r[end:-1:1,2],r[:,2], dims=1) .+ 1im .* cat(-r[end:-1:1,3],r[:,3], dims=1)
gt_h  = r[:,2] .+ 1im .* r[:,3]
#partial_sums = reverse([2*real(sum(gt[i:end-i+1]))/mP.β + 1 for i in 1:floor(Int,(length(gt)/2-100))]);
