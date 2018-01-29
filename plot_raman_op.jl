#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Utils: interactive
import NaCsCalc: Trap
using NaCsPlot
using Cubature
using PyPlot

matplotlib["rcParams"][:update](Dict("font.weight" => "normal"))

const m_Na = 23e-3 / 6.02e23
const m_Na = 23e-3 / 6.02e23
const η_rx = Trap.η(m_Na, 479e3, 2π / 589e-9) * √(2) * 0.96
const η_ry = Trap.η(m_Na, 492e3, 2π / 589e-9) * √(2) * 0.933
const η_az = Trap.η(m_Na, 85.7e3, 2π / 589e-9) * 0.67
const η_az_op = Trap.η(m_Na, 85.7e3, 2π / 589e-9)

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

function op_heating_az(n1, n2, η::T, isσ) where T
    hcubature(x->op_coupling(n1, n2, η, T(x[1]), T(x[2]), isσ, T(0), (T(1), T(0))),
              [0.0, 0.0], [1.0, π], abstol=1e-6)[1]
end

function op_heating_r(n1, n2, η::T, isσ) where T
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

const coupling_az = (op_heating_all(op_heating_az, 100, 100, Float32(η_az_op), false) *
                     op_heating_all(op_heating_az, 100, 100, Float32(η_az_op), true))
# const coupling_rx = (op_heating_all(op_heating_r, 30, 30, Float32(η_rx), false) *
#                      op_heating_all(op_heating_r, 30, 30, Float32(η_rx), true))
# const coupling_ry = (op_heating_all(op_heating_r, 30, 30, Float32(η_ry), false) *
#                      op_heating_all(op_heating_r, 30, 30, Float32(η_ry), true))

function plot_op_sidebands(n1s, n2s, coupling)
    for i in 1:length(n2s)
        n2 = n2s[i]
        plot(n1s, coupling[n2 + 1, n1s .+ 1], ".-")
    end
end

const nmax = 95
const nmax_r = 15

function plot_sidebands(ns, Δns, η)
    for Δn in Δns
        plot(ns, abs.(Trap.sideband.(ns, ns .+ Δn, η)), ".-", label="\$\\Delta n=$(Δn)\$")
    end
end

figure(figsize=[1.5, 1.1] * 4.8)

ax1 = subplot(211)
plot_op_sidebands(0:nmax, [0, 1, 2, 6, 14, 27, 45, 70], coupling_az)
text(0.5, 0.8, "\$n_{init}\\!\\!=\\!\\!0\$", color="C0")
text(1, 0.63, "\$n_{init}\\!\\!=\\!\\!1\$", color="C1")
text(2, 0.47, "\$n_{init}\\!\\!=\\!\\!2\$", color="C2")
text(3, 0.33, "\$n_{init}\\!\\!=\\!\\!6\$", color="C3")
text(7, 0.22, "\$n_{init}\\!\\!=\\!\\!14\$", color="C4")
text(17, 0.13, "\$n_{init}\\!\\!=\\!\\!27\$", color="C5")
text(38, 0.13, "\$n_{init}\\!\\!=\\!\\!45\$", color="C6")
text(60, 0.11, "\$n_{init}\\!\\!=\\!\\!70\$", color="C7")
text(85, 0.78, "(A)")
grid()
ylim([0, 0.9])
xlim([0, nmax])
ylabel("Probability")
setp(ax1[:get_xticklabels](), visible=false)
ax1[:tick_params](axis="x", length=0)
ax1[:get_yaxis]()[:set_label_coords](-0.105, 0.5)

ax2 = subplot(212)
subplots_adjust(hspace=0.05)
plot_sidebands(0:nmax, -1:-1:-5, η_az)
xlim([0, nmax])
ylim([0, 0.75])
grid()
ylabel("\$|\\langle n |e^{i\\Delta\\vec k\\cdot\\vec z}| n + \\Delta n \\rangle|\$")
text(0, 0.61, "\$\\Delta n\\!\\!=\\!\\!\\!-\\!1\$", color="C0")
text(12, 0.51, "\$\\Delta n\\!\\!=\\!\\!\\!-\\!2\$", color="C1")
text(29, 0.46, "\$\\Delta n\\!\\!=\\!\\!\\!-\\!3\$", color="C2")
text(47, 0.42, "\$\\Delta n\\!\\!=\\!\\!\\!-\\!4\$", color="C3")
text(74, 0.40, "\$\\Delta n\\!\\!=\\!\\!\\!-\\!5\$", color="C4")
text(85, 0.632, "(B)")
ax2[:get_yaxis]()[:set_label_coords](-0.105, 0.5)
xlabel("Motional state \$n\$")

NaCsPlot.maybe_save(joinpath(@__DIR__, "fig2_raman_op"))

NaCsPlot.maybe_show()
