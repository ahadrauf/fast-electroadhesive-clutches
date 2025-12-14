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
name_clarifier = "_release_time_vs_Fpreload"
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
file_name = "20241208_14_37_30_timing_vs_Fpreload_plotdata"

data = np.load(file_loc + file_name + ".npy", allow_pickle=True)
all_V, all_x1, all_y1, all_yerr1, all_x2, all_y2, all_yerr2, all_load_cell_engage_times, all_load_cell_disengage_times, all_load_cell_disengage_forces, all_load_cell_preload_forces, all_load_cell_initial_forces, all_load_cell_max_forces, all_load_cell_disengage_times_10pct = data

force_buckets = [(-1, 0.25), (0.25, 0.5), (0.5, 0.77), (0.77, 1.1), (1.1, 1.5)]
engagment_times_by_bucket = [[[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []]]
release_times_by_bucket = [[[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []]]
release_times_10pct_by_bucket = [[[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []]]
forces_by_bucket = [[[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []]]
for i, V in enumerate(all_V):
    all_load_cell_engage_times[V] = np.array(all_load_cell_engage_times[V])
    all_load_cell_disengage_forces[V] = np.array(all_load_cell_disengage_forces[V])
    all_load_cell_disengage_times[V] = np.array(all_load_cell_disengage_times[V])
    all_load_cell_initial_forces[V] = np.array(all_load_cell_initial_forces[V])
    all_load_cell_preload_forces[V] = np.array(all_load_cell_preload_forces[V])
    all_load_cell_disengage_times_10pct[V] = np.array(all_load_cell_disengage_times_10pct[V])

    for curr_te, curr_f in zip(all_load_cell_engage_times[V], all_load_cell_preload_forces[V]):
        for idx, bucket in enumerate(force_buckets):
            if bucket[0] <= curr_f <= bucket[1]:
                engagment_times_by_bucket[idx][i].append(curr_te)
                forces_by_bucket[idx][i].append(curr_f)

    for curr_tr, curr_f in zip(all_load_cell_disengage_times[V], all_load_cell_preload_forces[V]):
        for idx, bucket in enumerate(force_buckets):
            if bucket[0] <= curr_f <= bucket[1]:
                release_times_by_bucket[idx][i].append(curr_tr)

# Plot data
colors = ["#1964B0", "#F1932D", "#4DB264", "#DB060B"]  # "#F7F057" = yellow
fig, axs_all = plt.subplots(2, 2, layout='constrained', sharex='all', sharey='all', figsize=(6, 4))
axs = [axs_all[0][0], axs_all[0][1], axs_all[1][0], axs_all[1][1]]
print("Average force for bucket 0", np.mean([xi[0] for xi in all_x2]))

slopes_avg_data = []
slopes_raw_data = []
slopes_sim = []
min_x = np.min([np.mean(forces_by_bucket[0][i]) for i in range(len(all_V))])
max_x = np.max([np.mean(forces_by_bucket[-1][i]) for i in range(len(all_V))])
print("Min X, Max X", min_x, max_x)
for i, V in enumerate(all_V):
    x = [np.mean(fi[i]) for fi in forces_by_bucket]
    y = [np.mean(ti[i]) for ti in release_times_by_bucket]
    yerr = [np.std(ti[i]) for ti in release_times_by_bucket]

    slope, intercept, r, p, se = linregress(x, y)
    slopes_avg_data.append(slope)
    print("Linear fit data:", slope, intercept, r * r)
    slope, intercept, r, p, se = linregress(np.hstack([fi[i] for fi in forces_by_bucket]),
                                            np.hstack([ti[i] for ti in release_times_by_bucket]))
    slopes_raw_data.append(slope)
    print("Linear fit raw data:", slope, intercept, r * r)

    axs[i].errorbar(x, y, yerr=yerr,  # xerr=[np.std(fi[i]) for fi in forces_by_bucket],
                    capsize=5, ecolor='k', elinewidth=2, capthick=2, c='k')
    axs[i].scatter(x, y, c='k', s=50)

    loadcell_start_ratio_orig = np.mean(np.divide(all_load_cell_disengage_forces[V], all_load_cell_max_forces[V]))
    loadcell_start_ratio = np.mean(np.divide(all_load_cell_disengage_forces[V] - all_load_cell_initial_forces[V], all_load_cell_max_forces[V] - all_load_cell_initial_forces[V]))
    print("Average load cell start ratio:", loadcell_start_ratio, loadcell_start_ratio_orig)

    if True:
        Fpreload_range = np.linspace(min_x, max_x, 25)
        Fpreload_sim = []
        trelease_sim = []
        for Fpreload in Fpreload_range:
            L_dielectric = 55.5e-3
            w_pin = 2e-3
            depth_pin = 2e-3
            k_dielectric = 54.2
            t_dielectric = 24e-6
            frequency = 1000
            T, Y, dY, X, dX, events, t_release_10pct, _, _, _, _, _, _ = sim_release_time(L_dielectric=L_dielectric,
                                                                                          w_pin=w_pin, depth_pin=depth_pin,
                                                                                          V_max=V, t_dielectric=t_dielectric,
                                                                                          frequency=frequency, Fext=Fpreload,
                                                                                          loadcell_start_ratio=loadcell_start_ratio)
            if type(events) == np.float64:
                Fpreload_sim.append(Fpreload)
                trelease_sim.append(events)
        axs[i].plot(Fpreload_sim, trelease_sim, ls='--', c='tab:red', lw=2)
        slope, intercept, r, p, se = linregress(Fpreload_sim, trelease_sim)
        print("Linear fit model:", slope, intercept, r * r)
        slopes_sim.append(slope)

    bbox = dict(alpha=0.8, boxstyle='round', fc='white', ec='0.8')
    axs[i].text(0.05, 0.93, "{} V".format(int(V)), bbox=bbox, transform=axs[i].transAxes, va='top', horizontalalignment='left')
    axs[i].grid(True)
    axs[i].set_xticks([0, 0.5, 1, 1.5])
    axs[i].yaxis.set_tick_params(labelleft=True)
    axs[i].xaxis.set_tick_params(labelbottom=True)

fig.supxlabel("Preload Normal Force (N)", fontsize=16)
fig.supylabel("Release Time (ms)", fontsize=16)
fig.text(0.00, .98, "(d)", transform=fig.transFigure, horizontalalignment='left', verticalalignment='top')

# plt.savefig(save_folder + timestamp + ".png", dpi=300)
# # # plt.savefig("figures/" + timestamp + ".svg")
# plt.savefig(save_folder + timestamp + ".pdf")

plt.show()
