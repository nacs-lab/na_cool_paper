#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))
import NaCsCalc: Trap

# Matrix elements
const m_Na = 23e-3 / 6.02e23
const νs = (69e3, 430e3, 590e3)
const ηs_op = Trap.η.(m_Na, νs, 2π / 589e-9)
const ηs_ra = ηs_op .* (0.67, √2, √2)

function calc_nbar(TμK)
    ħ = 1.0545718e-34
    k_B = 1.38064852e-23
    Hz_per_μK = k_B * 1e-6 / ħ / 2π
    THz = TμK * Hz_per_μK
    return THz ./ νs
end

function calc_η_eff(TμK)
    return sqrt.(2 .* calc_nbar(TμK) .+ 1) .* ηs_op
end

@show ηs_op
@show ηs_ra
@show calc_nbar(70)
@show calc_η_eff(70)
# @show calc_η_eff(75)
