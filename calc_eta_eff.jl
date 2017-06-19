#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))
import NaCsCalc: Trap

# Matrix elements
const m_Na = 23e-3 / 6.02e23
const ν_a1 = 68.8e3
const ν_r2 = 430e3
const ν_r3 = 589.5e3
const η_a1_op = Trap.η(m_Na, ν_a1, 2π / 589e-9)
const η_r2_op = Trap.η(m_Na, ν_r2, 2π / 589e-9)
const η_r3_op = Trap.η(m_Na, ν_r3, 2π / 589e-9)
const η_a1 = η_a1_op * 0.67
const η_r2 = η_r2_op * √(2)
const η_r3 = η_r3_op * √(2)

function calc_η_eff(TμK)
    ħ = 1.0545718e-34
    k_B = 1.38064852e-23
    Hz_per_μK = k_B * 1e-6 / ħ / 2π
    THz = TμK * Hz_per_μK
    nbars = THz ./ (ν_a1, ν_r2, ν_r3)
    return sqrt.(2 .* nbars .+ 1) .* (η_a1_op, η_r2_op, η_r3_op)
end

@show calc_η_eff(70)
