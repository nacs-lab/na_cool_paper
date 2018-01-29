#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "../../experiments/misc/data/data_20170402_205344.csv")
const iname_b = joinpath(@__DIR__, "../../experiments/misc/data/data_20170404_133229.csv")
const iname_c = joinpath(@__DIR__, "../../experiments/misc/data/data_20170409_005850.csv")
const iname_d = joinpath(@__DIR__, "../../experiments/misc/data/data_20170409_082523.csv")

const data_a = NaCsData.load_count_csv(iname_a)
const data_b = NaCsData.load_count_csv(iname_b)
const data_c = NaCsData.load_count_csv(iname_c)
const data_d = NaCsData.load_count_csv(iname_d)

const spec_a = OrderedDict(
    # With cooling +-1
    :cool_pm1=>((linspace(-18.985, -18.785, 11), linspace(-18.15, -17.95, 11)),
                (linspace(-18.945, -19.145, 11), linspace(-18.00, -17.75, 11)),
                (linspace(-18.545, -18.585, 11), linspace(-18.415, -18.455, 11))),
    # Without cooling +-1, -2
    :nocool_pm12=>((linspace(-18.985, -18.785, 11), linspace(-18.15, -17.95, 11),
                    linspace(-17.53, -17.745, 11)),
                   (linspace(-18.945, -19.145, 11), linspace(-18.00, -17.75, 11),
                    linspace(-17.14, -17.44, 11)),
                   (linspace(-18.545, -18.585, 11), linspace(-18.415, -18.455, 11),
                    linspace(-18.35, -18.39, 11))),
    # Without cooling carrier
    :nocool_0=>(linspace(-18.335, -18.58, 11),
                linspace(-18.31, -18.61, 11),
                linspace(-18.47, -18.53, 11),),
    # Without cooling axial high orders
    :nocool_a8=>linspace(-18.35, -17.90, 91),
    # Without repeating
    :norepeat=>((linspace(-18.985, -18.785, 11), linspace(-18.15, -17.95, 11)),
                (linspace(-18.945, -19.145, 11), linspace(-18.00, -17.75, 11)),
                (linspace(-18.545, -18.585, 11), linspace(-18.415, -18.455, 11))),
    # With waiting
    :wait=>((linspace(-18.985, -18.785, 11), linspace(-18.15, -17.95, 11)),
            (linspace(-18.945, -19.145, 11), linspace(-18.00, -17.75, 11)),
            (linspace(-18.545, -18.585, 11), linspace(-18.415, -18.455, 11)))
)
const spec_b = (linspace(-17.53, -17.745, 11),
                linspace(-17.14, -17.44, 11),
                linspace(-18.40, -17.90, 51),
                linspace(-18.40, -18.35, 11),)
const spec_c = ((linspace(-18.985, -18.785, 11), linspace(-18.15, -17.95, 11),
                 linspace(-17.53, -17.745, 11)),
                (linspace(-18.945, -19.145, 11), linspace(-18.00, -17.75, 11),
                 linspace(-17.14, -17.44, 11)),
                (linspace(-18.54, -18.58, 11), linspace(-18.400, -18.440, 11)),
                # Axial -2 ~ -8
                linspace(-18.40, -17.90, 51) # 4
                )
const spec_d = (
    (linspace(-18.985, -18.785, 11), linspace(-18.15, -17.95, 11),
     linspace(-17.53, -17.745, 11)),
    (linspace(-18.945, -19.145, 11), linspace(-18.00, -17.75, 11),
     linspace(-17.14, -17.44, 11)),
    (linspace(-18.54, -18.58, 11), linspace(-18.400, -18.440, 11))
)

const split_a = NaCsData.split_data(data_a, spec_a)
const split_b = NaCsData.split_data(data_b, spec_b)
const split_c = NaCsData.split_data(data_c, spec_c)
const split_d = NaCsData.split_data(data_d, spec_d)

const prefix = joinpath(@__DIR__, "spectrum")

to_sideband(f) = (i, v)->(v - f) * 1000

data_nocool_r2 = NaCsData.map_params(to_sideband(-18.4625), split_a[:nocool_pm12][1])
data_nocool_r2_0 = NaCsData.map_params(to_sideband(-18.4625), split_a[:nocool_0][1])
data_cool_r2 =  NaCsData.map_params(to_sideband(-18.4575), split_d[1])

data_nocool_r3 = NaCsData.map_params(to_sideband(-18.485), split_a[:nocool_pm12][2])
data_nocool_r3_0 = NaCsData.map_params(to_sideband(-18.485), split_a[:nocool_0][2])
data_cool_r3 =  NaCsData.map_params(to_sideband(-18.480), split_d[2])

data_nocool_a1 = NaCsData.map_params(to_sideband(-18.4965), split_a[:nocool_pm12][3])
data_nocool_a1_0 = NaCsData.map_params(to_sideband(-18.4965), split_a[:nocool_0][3])
data_nocool_a1_hi = [NaCsData.map_params(to_sideband(-18.5015), split_a[:nocool_a8]);
                     NaCsData.map_params(to_sideband(-18.4965), split_b[4])]
data_cool_a1 =  NaCsData.map_params(to_sideband(-18.488), split_d[3])
data_cool_a1_hi =  NaCsData.map_params(to_sideband(-18.488), split_c[4])

fig = figure(figsize=[1.6, 1] * 4.8)

# Without cooling
# Radial 2
NaCsPlot.plot_survival_data(data_nocool_r2[1], fmt="C3^-", label="\$x\$ initial")
NaCsPlot.plot_survival_data(data_nocool_r2[2], fmt="C3^-")
NaCsPlot.plot_survival_data(data_nocool_r2[3], fmt="C3^-")
# Radial 3
NaCsPlot.plot_survival_data(data_nocool_r3[1], fmt="C3o--", label="\$y\$ initial")
NaCsPlot.plot_survival_data(data_nocool_r3[2], fmt="C3o--")
NaCsPlot.plot_survival_data(data_nocool_r3[3], fmt="C3o--")
text(-660, 0.9, "(A)")

# With cooling
# Radial 2
NaCsPlot.plot_survival_data(data_cool_r2[1], fmt="C0v-", label="\$x\$ cooled")
NaCsPlot.plot_survival_data(data_cool_r2[2], fmt="C0v-")
NaCsPlot.plot_survival_data(data_cool_r2[3], fmt="C0v-")
# Radial 3
NaCsPlot.plot_survival_data(data_cool_r3[1], fmt="C0s--", label="\$y\$ cooled")
NaCsPlot.plot_survival_data(data_cool_r3[2], fmt="C0s--")
NaCsPlot.plot_survival_data(data_cool_r3[3], fmt="C0s--")
grid()
ylim([0, 1])
xlim([-700, 1500])
legend()
xlabel("\$\\delta\$, Detuning from carrier (kHz)")
ylabel("F=1 population")

NaCsPlot.maybe_save("$(prefix)_r")

fig = figure(figsize=[1.6, 1] * 4.8)
# Without cooling
NaCsPlot.plot_survival_data(data_nocool_a1[1], fmt="C3o-", label="\$z\$ initial")
NaCsPlot.plot_survival_data(data_nocool_a1[2], fmt="C3o-")
NaCsPlot.plot_survival_data(data_nocool_a1_0, fmt="C3o-")
NaCsPlot.plot_survival_data(data_nocool_a1_hi, fmt="C3o-")

# With cooling
NaCsPlot.plot_survival_data(data_cool_a1[1], fmt="C0s-", label="\$z\$ cooled")
NaCsPlot.plot_survival_data(data_cool_a1[2], fmt="C0s-")
NaCsPlot.plot_survival_data(data_cool_a1_hi, fmt="C0s-")
grid()
ylim([0, 0.6])
xlim([-100, 620])
legend()
text(-86, 0.54, "(A)")
xlabel("\$\\delta\$, Detuning from carrier (kHz)")
ylabel("F=1 population")
NaCsPlot.maybe_save("$(prefix)_a1")

NaCsPlot.maybe_show()
