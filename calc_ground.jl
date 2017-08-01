#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))
import NaCsCalc.Format: Unc

const a1m1 = [0.0410, 0.0182]
const a1p1 = [0.476, 0.411]

const r2m1 = [0.0656, 0.0366]
const r2p1 = [0.656, 0.593]

const r3m1 = [0.0310, 0.0121]
const r3p1 = [0.459, 0.395]

function naive_ground(m1, p1)
    avg_m1 = (m1[1] + m1[2]) / 2
    std_m1 = abs(m1[1] - m1[2]) / 2
    avg_p1 = (p1[1] + p1[2]) / 2
    std_p1 = abs(p1[1] - p1[2]) / 2

    avg_r = avg_m1 / avg_p1
    std_r = avg_r * sqrt((std_p1 / avg_p1)^2 + (std_m1 / avg_m1)^2)

    return Unc(1 - avg_r, std_r)
end

@show naive_ground(a1m1, a1p1)
@show naive_ground(r2m1, r2p1)
@show naive_ground(r3m1, r3p1)
