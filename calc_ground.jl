#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))
import NaCsCalc.Format: Unc

const azm1 = [0.0410, 0.0182]
const azp1 = [0.476, 0.411]

const rxm1 = [0.0656, 0.0366]
const rxp1 = [0.656, 0.593]

const rym1 = [0.0310, 0.0121]
const ryp1 = [0.459, 0.395]

const azm10 = [0.419, 0.317]
const azp10 = [0.497, 0.396]

const rxm10 = [0.599, 0.500]
const rxp10 = [0.820, 0.740]

const rym10 = [0.455, 0.355]
const ryp10 = [0.664, 0.574]

function naive_ground(m1, p1)
    avg_m1 = (m1[1] + m1[2]) / 2
    std_m1 = abs(m1[1] - m1[2]) / 2
    avg_p1 = (p1[1] + p1[2]) / 2
    std_p1 = abs(p1[1] - p1[2]) / 2

    avg_r = avg_m1 / avg_p1
    std_r = avg_r * sqrt((std_p1 / avg_p1)^2 + (std_m1 / avg_m1)^2)

    avg_g = 1 - avg_r
    std_g = std_r

    avg_n = 1 / avg_g - 1
    std_n = std_g / avg_g^2

    return Unc(avg_g, std_g), Unc(avg_n, std_n)
end

function show_ground(m1, p1, name)
    println(name)
    g, n = naive_ground(m1, p1)
    println("  P_0: ", g)
    println("  nbar: ", n)
end

show_ground(rxm1, rxp1, "Radial (x) cooled")
show_ground(rym1, ryp1, "Radial (y) cooled")
show_ground(azm1, azp1, "Axial (z) cooled")

show_ground(rxm10, rxp10, "Radial (x) initial")
show_ground(rym10, ryp10, "Radial (y) initial")
show_ground(azm10, azp10, "Axial (z) initial")
