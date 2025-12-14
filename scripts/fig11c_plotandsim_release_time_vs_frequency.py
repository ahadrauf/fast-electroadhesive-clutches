import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit, root_scalar
from scipy.integrate import solve_ivp
from scipy.special import kv
from scipy.stats import linregress
from datetime import datetime
from sim_release_time import *

now = datetime.now()
name_clarifier = "_release_time_vs_frequency_plotandsim"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
save_folder = "../figures/"

plt.rcParams["font.family"] = "Arial"
plt.rc('font', size=16)  # controls default text size
plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=16)  # fontsize of the y tick labels
plt.rc('legend', fontsize=16)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

file_loc = "../data/"
file_name = "20241208_14_37_18_timing_vs_frequency_plotandsim_plotdata"

data = np.load(file_loc + file_name + ".npy", allow_pickle=True)
all_V, all_x1, all_y1, all_yerr1, all_x2, all_y2, all_yerr2, all_load_cell_engage_times, all_load_cell_disengage_times, all_load_cell_disengage_forces, all_load_cell_preload_forces, all_load_cell_initial_forces, all_load_cell_max_forces, all_load_cell_disengage_times_10pct = data
all_freq = sorted(all_load_cell_preload_forces[150].keys())

# Plot data
colors = ["#1964B0", "#F1932D", "#4DB264", "#DB060B"]  # "#F7F057" = yellow
fig, axs_all = plt.subplots(2, 2, layout='constrained', sharex='all', sharey='all',
                            figsize=(6, 4))
axs = [axs_all[0][0], axs_all[0][1], axs_all[1][0], axs_all[1][1]]
slopes_avg_data = []
slopes_raw_data = []
slopes_sim = []
tau_sim = []
ratios_1Hz = []
ratios_10kHz = []
ratios = []
for i, V in enumerate(all_V):
    x = all_freq
    y = [np.mean(all_load_cell_disengage_times[V][freq]) for freq in all_freq]
    yerr = [np.std(all_load_cell_disengage_times[V][freq]) for freq in all_freq]

    slope, intercept, r, p, se = linregress(np.log10(x), y)
    slopes_avg_data.append(slope)
    print("Linear fit data:", slope, intercept, r*r)
    slope, intercept, r, p, se = linregress(np.hstack([[np.log10(freq)]*len(all_load_cell_disengage_times[V][freq]) for freq in all_freq]),
                                            np.hstack([all_load_cell_disengage_times[V][freq] for freq in all_freq]))
    slopes_raw_data.append(slope)
    print("Linear fit raw data:", slope, intercept, r*r)
    print(V, x)
    print(V, y)
    ratios.append(y[1]/y[4])  # ratio between 1 Hz and 1 kHz
    print(y[1]/y[4], y[1], y[4])

    errorbar = axs[i].errorbar(x, y, yerr=yerr,  # xerr=[np.std(fi[i]) for fi in forces_by_bucket],
                               capsize=5, ecolor='k', elinewidth=2, capthick=2, c='k')  # alpha=0.4
    scatter = axs[i].scatter(x, y, c='k', s=50)

    frequency_range = np.power(10, np.linspace(-1, 4, 10))

    freq_sim = []
    tengage_sim = []
    all_preloads = np.hstack([all_load_cell_preload_forces[V][f] for f in all_freq])
    all_initial_forces = np.hstack([all_load_cell_initial_forces[V][f] for f in all_freq])
    all_disengage_forces = np.hstack([all_load_cell_disengage_forces[V][f] for f in all_freq])
    all_max_forces = np.hstack([all_load_cell_max_forces[V][f] for f in all_freq])
    print("Preloads for V =", V, np.mean(all_preloads), all_preloads)
    Fpreload = np.mean(all_preloads)

    loadcell_start_ratio = np.mean(np.divide(all_disengage_forces - all_initial_forces, all_max_forces - all_initial_forces))
    Fpreloads_by_frequecy = [np.mean(all_load_cell_preload_forces[V][f]) for f in all_freq]
    loadcell_start_ratios_by_frequency = [np.mean(np.divide(all_load_cell_disengage_forces[V][f], all_load_cell_max_forces[V][f])) for f in all_freq]

    print("Ratio between 10 Hz and 1 kHz", y[-4]/y[1], y[1], y[-4])
    print("Ratio between 10 kHz and 1 kHz", y[-4]/y[-1], y[2], y[-4])
    ratios_1Hz.append(y[-3]/y[1])
    ratios_10kHz.append(y[-3]/y[-1])

    if True:
        for frequency in frequency_range:
            L_dielectric = 55.5e-3
            w_pin = 2e-3
            depth_pin = 2e-3
            k_dielectric = 54.2
            t_dielectric = 24e-6

            # Screenshot 2024-12-08 190214 = both enabled, loadcell only made no difference, Screenshot 2024-12-08 191525 = preload only
            Fpreload = np.interp(np.log10(frequency), np.log10(all_freq), Fpreloads_by_frequecy)
            T, Y, dY, X, dX, events, t_release_10pct, _, _, _, _, _, _ = sim_release_time(L_dielectric=L_dielectric,
                                                                                          w_pin=w_pin, depth_pin=depth_pin,
                                                                                          V_max=V, t_dielectric=t_dielectric,
                                                                                          frequency=frequency, Fext=Fpreload,
                                                                                          loadcell_start_ratio=loadcell_start_ratio)
            if type(events) == np.float64:
                freq_sim.append(frequency)
                tengage_sim.append(events)
        axs[i].semilogx(freq_sim, tengage_sim, ls='--', c='tab:red', lw=2)
        slope, intercept, r, p, se = linregress(np.log10(freq_sim), tengage_sim)
        print("Linear fit model:", slope, intercept, r*r)
        slopes_sim.append(slope)

    bbox = dict(alpha=0.8, boxstyle='round', fc='white', ec='0.8')
    axs[i].text(0.95, 0.93, "{} V".format(int(V)), bbox=bbox, transform=axs[i].transAxes, horizontalalignment='right', verticalalignment='top')
    axs[i].grid(True)
    axs[i].tick_params(labelleft=True)
    axs[i].set_xscale('log')
    axs[i].set_xticks([0.1, 1, 10, 100, 1000, 10000], ["DC", r"$10^0$", r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$"])
    axs[i].xaxis.set_tick_params(labelbottom=True)

fig.supxlabel("Drive Frequency (Hz)", fontsize=16)
fig.supylabel("Release Time (ms)", fontsize=16)
fig.text(0.00, .98, "(c)", transform=fig.transFigure, horizontalalignment='left', verticalalignment='top')

print("All Slopes Avg", np.mean(slopes_avg_data), slopes_avg_data)
print("All Slopes Raw", np.mean(slopes_raw_data), slopes_raw_data)
print("All Slopes Sim", np.mean(slopes_sim), slopes_sim)
print("Ratios 1 Hz", np.mean(ratios_1Hz), ratios_1Hz)
print("Ratios 10 kHz", np.mean(ratios_10kHz), ratios_10kHz)
print("Ratios between 1 Hz and 1 kHz", np.mean(ratios), ratios)

# file_name = "20241214_05_46_15_timing_vs_frequency_plotandsim_tr10pct_plotdata"
# data = np.load(save_folder + file_name + ".npy", allow_pickle=True)
# all_load_cell_disengage_times_10pct = data[-1]
#
# fig2, axs2_all = plt.subplots(2, 2, layout='constrained', sharex='all', sharey='all', figsize=(6, 4))
# axs2 = axs2_all.flat
# for i, V in enumerate(all_V):
#     x = all_freq
#     y = [np.mean(all_load_cell_disengage_times_10pct[V][freq]) for freq in all_freq]
#     yerr = [np.std(all_load_cell_disengage_times_10pct[V][freq]) for freq in all_freq]
#     errorbar = axs2[i].errorbar(x, y, yerr=yerr, capsize=5, ecolor='k', elinewidth=2, capthick=2, c='k', ls='--')  # alpha=0.4
#     scatter = axs2[i].scatter(x, y, c='k', s=50, marker='x')
#
#     x = all_freq
#     y = [np.mean(all_load_cell_disengage_times[V][freq]) for freq in all_freq]
#     yerr = [np.std(all_load_cell_disengage_times[V][freq]) for freq in all_freq]
#     errorbar4 = axs2[i].errorbar(x, y, yerr=yerr, capsize=5, ecolor='k', elinewidth=2, capthick=2, color='k', lw=2)
#     scatter4 = axs2[i].scatter(x, y, c='k', s=50, marker='o')
#
#     axs2[i].text(0.95, 0.93, "{} V".format(int(V)), bbox=bbox, transform=axs2[i].transAxes, horizontalalignment='right', verticalalignment='top')
#     axs2[i].grid(True)
#     axs2[i].tick_params(labelleft=True)
#     axs2[i].set_xscale('log')
#     axs2[i].set_xticks([0.1, 1, 10, 100, 1000, 10000], ["DC", r"$10^0$", r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$"])
#     axs2[i].xaxis.set_tick_params(labelbottom=True)
#     # if i == 1:
#     #     # axs2[i].legend([(scatter4, errorbar4), (scatter, errorbar)], [r"$t_{r,10\%}$", r"$t_{r,90\%}$"])
#     #     axs2[i].legend([(scatter4, errorbar4), (scatter, errorbar)], ["90% Fall Time", "10% Fall Time"])
# fig2.supxlabel("Drive Frequency (Hz)", fontsize=16)
# # fig2.supylabel("Release Time (ms)", fontsize=16)
# fig2.supylabel("Load Cell Fall Times (ms)", fontsize=16)
# fig2.text(0.00, .98, "(c)", transform=fig2.transFigure, horizontalalignment='left', verticalalignment='top')

# plt.savefig(save_folder + timestamp + ".png", dpi=300)
# # # plt.savefig("figures/" + timestamp + ".svg")
# plt.savefig(save_folder + timestamp + ".pdf")
# fig2.savefig(save_folder + timestamp + ".png", dpi=300)
# # # plt.savefig("figures/" + timestamp + ".svg")
# fig2.savefig(save_folder + timestamp + ".pdf")

plt.show()
