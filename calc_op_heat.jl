#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc: Trap
using Cubature

const m_Na = 23e-3 / 6.02e23
const η_ax = Trap.η(m_Na, 67e3, 2π / 589e-9) / √(2)
const η_ax_op = Trap.η(m_Na, 67e3, 2π / 589e-9)

"""
    emission(v, isσ::Bool) -> (sinθ, cosθ)

Given an uniform distribution in `[0, 1]` of `v` this function will return a distribution
of `cosθ` that matches a dipole emission pattern (the type of the transition is determined by
`isσ`).
"""
function emission_angle(v::T, isσ::Bool) where T<:AbstractFloat
    # Returns `cosθ`. The caller should be ready to handle `|cosθ| > 1`
    # The PDFs of the `θ` distribution are
    # `3 / 4 * (1 - cos²θ) * sinθ` for π light and
    # `3 / 8 * (1 + cos²θ) * sinθ` for σ± light
    # The corresponding CDFs are
    # `1 / 4 * (3cosθ - cos³θ) + 1 / 2` for π light and
    # `1 / 8 * (3cosθ + cos³θ) + 1 / 2` for σ± light
    # For a random number `v` we need to solve
    # `3cosθ - cos³θ = 4v - 2` for π light and
    # `3cosθ + cos³θ = 8v - 4` for σ± light
    # The (real) solution is
    # `cosθ = x + 1 / x` where `x = ∛(2√(v² - v) - 2v + 1)` for π light and
    # `cosθ = x - 1 / x` where `x = ∛(√(16v² - 16v + 5) + 4v - 2)` for σ± light
    # Note that the `x` for π polarization is complex
    if isσ
        y = muladd(v, 4, -2)
        x = @fastmath cbrt(sqrt(muladd(y, y, 1)) + y)
        return x - 1 / x
    else
        θ′ = @fastmath acos(muladd(T(-2), v, T(1)))
        cosθ = @fastmath cos(muladd(θ′, T(1 / 3), - T(2π / 3))) * 2
        return cosθ
    end
end

"""
    op_coupling(n1, n2, η::T, v::T, φ::T, isσ::Bool, cosθ_dri::T, sincosθ_quan::T)

Given an uniform distribution in `[0, 1]` of `v` and an uniform distribution in `[0, π]` of φ
this function will return the distribution of matrix element of matrix element between `n1`
and `n2` motional states after one photon scattering.

The polarization of the emission is determined by `isσ`.
The overlap between the axis and the drive beam direction is `cosθ_dri`.
The `sincos` of the angle between the axis and the quantization axis is `sincosθ_quan`.
"""
function op_coupling(n1, n2, η::T, v::T, φ::T, isσ::Bool, cosθ_dri::T,
                     sincosθ_quan::Tuple{T,T}) where T
    cosθ::T = emission_angle(v, isσ)
    sinθ::T = if abs(cosθ) >= 1
        T(0)
    else
        @fastmath sqrt(1 - cosθ^2)
    end
    sinθ_quan, cosθ_quan = sincosθ_quan
    η_eff = @fastmath η * (cosθ * cosθ_quan + sinθ * sinθ_quan * cos(φ) + cosθ_dri)
    return Trap.sideband(n1, n2, abs(η_eff))^2 / π
end

function op_heating_ax(n1, n2, η::T, isσ) where T
    hcubature(x->op_coupling(n1, n2, η, T(x[1]), T(x[2]), isσ, T(0), (T(1), T(0))),
              [0.0, 0.0], [1.0, π], abstol=1e-6)[1]
end

function op_heating_ra(n1, n2, η::T, isσ) where T
    hcubature(x->op_coupling(n1, n2, η, T(x[1]), T(x[2]), isσ,
                             sqrt(T(0.5)), (sqrt(T(0.5)), sqrt(T(0.5)))),
              [0.0, 0.0], [1.0, π], abstol=1e-6)[1]
end

function op_heating_all(cb::Function, sz1, sz2, η::T, isσ) where T
    res = Matrix{T}(sz1, sz2)
    @inbounds for j in 1:sz2
        for i in 1:sz1
            fill_reverse = false
            if i < j && sz1 > sz2
                continue
            elseif i > j && sz1 <= sz2
                continue
            end
            v = cb(i - 1, j - 1, η, isσ)
            res[i, j] = v
            if i <= sz2 && j <= sz1
                res[j, i] = v
            end
        end
    end
    res
end

function p_heat(coupling, sz)
    T = eltype(coupling)
    res = Vector{T}(sz)
    @inbounds for i in 1:sz
        v = zero(T)
        for j in (i + 1):size(coupling, 1)
            v += coupling[j, i]
        end
        res[i] = v
    end
    return res
end

function op_unc(coupling, sz)
    T = eltype(coupling)
    res = Vector{T}(sz)
    @inbounds for i in 1:sz
        sump = zero(T)
        avgΔn = zero(T)
        @simd for j in 1:size(coupling, 1)
            p = coupling[j, i]
            sump += p
            avgΔn += abs(j - i) * p
        end
        res[i] = avgΔn / sump
    end
    return res
end

function prob_nbar(n, nbar)
    α = nbar / (nbar + 1)
    return α^n / (nbar + 1)
end

function calc_heating(cb::F, nbar, η, nσ) where F
    nmax = ceil(Int, nbar * 4)
    cpl_max = ceil(Int, nbar * 5)
    cplπ = op_heating_all(cb, cpl_max, cpl_max, Float32(η), false)
    cplσ = op_heating_all(cb, cpl_max, cpl_max, Float32(η), false)
    cpl = cplπ
    for _ in 1:nσ
        cpl *= cplσ
    end
    # heating = p_heat(cpl, nmax + 1)
    heating = op_unc(cpl, nmax + 1)
    probability = prob_nbar.(0:nmax, nbar)
    nbar_lo = floor(Int, nbar)
    w_hi = nbar - nbar_lo
    w_lo = 1 - w_hi
    heating_mean = heating[nbar_lo + 1] * w_lo + heating[nbar_lo + 2] * w_hi
    return sum(heating .* probability), heating[1], heating_mean
end

function prn_heating_ax(prefix, nbar, η, nσ)
    heat_avg, heat0, heat_nbar = calc_heating(op_heating_ax, nbar, η, nσ)
    println("$prefix:")
    println("  Avg heating: $heat_avg")
    println("  Heating for ground state: $heat0")
    println("  Heating for nbar: $heat_nbar")
end

function prn_heating_ra(prefix, nbar, η, nσ)
    heat_avg, heat0, heat_nbar = calc_heating(op_heating_ra, nbar, η, nσ)
    println("$prefix:")
    println("  Avg heating: $heat_avg")
    println("  Heating for ground state: $heat0")
    println("  Heating for nbar: $heat_nbar")
end

prn_heating_ax("Na", 20.9, η_ax_op, 1)
prn_heating_ax("Cs projection Ax", 6.0, 0.37, 2)
prn_heating_ra("Cs projection Ra", 6.0, 0.37, 2)
prn_heating_ax("Cs", 9.2, 0.32, 2)
prn_heating_ax("Rb", 8, 0.34, 2)
