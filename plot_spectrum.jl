#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using PyPlot
using DataStructures
matplotlib["rcParams"][:update](Dict("font.size" => 20,
                                     "font.weight" => "bold"))
matplotlib[:rc]("xtick", labelsize=15)
matplotlib[:rc]("ytick", labelsize=15)

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

function plot_data(data, scale=1; kws...)
    params, ratios, uncs = NaCsData.get_values(data)
    perm = sortperm(params)
    params = params[perm]
    ratios = ratios[perm, 2] .* scale
    uncs = uncs[perm, 2] .* scale
    errorbar(params, ratios, uncs; kws...)
end

function maybe_save(name)
    if !interactive()
        savefig("$name.png"; bbox_inches="tight", transparent=true)
        savefig("$name.svg", bbox_inches="tight", transparent=true)
        close()
    end
end

function maybe_show()
    if interactive()
        show()
    end
end

const prefix = joinpath(@__DIR__, "imgs/spectrum")

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

fig = figure(figsize=[1.5, 1] * 4.8)

# Without cooling
ax1 = subplot(211)
# Radial 2
plot_data(data_nocool_r2[1], 1 / 0.95, fmt="ro-", label="Axis 2")
plot_data(data_nocool_r2[2], 1 / 0.95, fmt="ro-")
plot_data(data_nocool_r2[3], 1 / 0.95, fmt="ro-")
# Radial 3
plot_data(data_nocool_r3[1], 1 / 0.95, fmt="ko-", label="Axis 3")
plot_data(data_nocool_r3[2], 1 / 0.95, fmt="ko-")
plot_data(data_nocool_r3[3], 1 / 0.95, fmt="ko-")
grid()
ylim([0, 1])
xlim([-700, 1400])
legend()
setp(ax1[:get_xticklabels](), visible=false)

# With cooling
ax2 = subplot(212)
subplots_adjust(hspace=0)
axvline(605, linewidth=4, color="k", ls="--", alpha=0.3, ymin=0.08)
axvline(1210, linewidth=4, color="k", ls="--", alpha=0.3, ymin=0.05)
axvline(430, linewidth=4, color="r", ls="--", alpha=0.3, ymin=0.1)
axvline(860, linewidth=4, color="r", ls="--", alpha=0.3, ymin=0.05)
# Radial 2
plot_data(data_cool_r2[1], 1 / 0.85, fmt="ro-")
plot_data(data_cool_r2[2], 1 / 0.85, fmt="ro-")
plot_data(data_cool_r2[3], 1 / 0.85, fmt="ro-")
# Radial 3
plot_data(data_cool_r3[1], 1 / 0.85, fmt="ko-")
plot_data(data_cool_r3[2], 1 / 0.85, fmt="ko-")
plot_data(data_cool_r3[3], 1 / 0.85, fmt="ko-")
grid()
ylim([0, 0.95])
xlim([-700, 1400])
fig[:text](0.5, 0.0, "Detuning from carrier (kHz)", ha="center", va="center")
fig[:text](0.0, 0.5, "Normalized survival", ha="center", va="center", rotation="vertical")

maybe_save("$(prefix)_r")

# figure(figsize=[2.5, 1] * 4.8)
# # Without cooling
# plot_data(data_nocool_a1[1], 1 / 0.95, fmt="ro-", label="Before")
# plot_data(data_nocool_a1[2], 1 / 0.95, fmt="ro-")
# plot_data(data_nocool_a1_0, 1 / 0.95, fmt="ro-")
# plot_data(data_nocool_a1_hi, 1 / 0.95, fmt="r^-", label="Before")
# # With cooling
# plot_data(data_cool_a1[1], 1 / 0.85, fmt="bo-", label="After")
# plot_data(data_cool_a1[2], 1 / 0.85, fmt="bo-")
# plot_data(data_cool_a1_hi, 1 / 0.85, fmt="b^-", label="After")
# grid()
# ylim([0, 0.8])
# xlim([-100, 620])
# title("Axis 1 (axial)")
# xlabel("Detuning from carrier (kHz)")
# ylabel("Normalized survival")
# legend()
# maybe_save("$(prefix)_a1")

maybe_show()
