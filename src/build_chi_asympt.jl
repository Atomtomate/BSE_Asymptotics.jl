"""

"""
function build_Fasympt!(Fasympt, Fcore, χ_sp, χ_ch, χ_pp, λ_ch, χ₀, U::Float64)
    # TODO: build border with λ
    # TODO: build corners with χ
    for ωi in axes(Fasympt,3)
        #TODO: ν, νp loops
        for νi in TODO_nu_list
            for νpi in TODO_nup_list
                #TODO: index to freq
            end
        end
    end
end

function calc_χλ!(χ, λ, F, χ₀, U::Float64, β::Float64)
    shape = size(F)                         # axis: ν. ν', ω
    F_tmp = similar(F, shape[1], shape[2])  # ν x ν'

    for ωi in axes(F,3)
        fill!(F_tmp, zero(eltype(F)))
        for qi in axes(χ₀,1)
            for νi in axes(F_tmp,1)
                F_tmp[νi,νi] = χ₀[qi,νi,ωi]
                for νpi in axes(F_tmp,2)
                    F_tmp[νi,νpi] -= χ₀[qi,νi,ωi]*F[νi,νpi,ωi]*χ₀[qi,νpi,ωi]/(β^2)
                end
            end
            χ[qi, ωi] = sum(F_tmp)/β
            for νk in axes(F, 3)
                λ[qi, νk, ωi] = sum(F_tmp[:,νk]) / (χ₀[qi, νk, ωi] * (1.0 + U * χ[qi, ωi]))
            end
        end
    end
    return χ, λ
end
