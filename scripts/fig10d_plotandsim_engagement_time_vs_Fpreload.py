import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit, root_scalar
from scipy.integrate import solve_ivp
from scipy.special import kv
from datetime import datetime
import csv
# from sim_engagement_time_20241121 import *
from sim_engagement_time import *

now = datetime.now()
name_clarifier = "_engagement_time_vs_Fpreload_plotandsim"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
save_folder = "C:/Users/ahadrauf/Desktop/Research/latex/electroadhesive_dynamics_paper/figures_test/"

plt.rcParams["font.family"] = "Arial"
# plt.rcParams["font.family"] = "serif"
# plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rc('font', size=16)  # controls default text size
plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=16)  # fontsize of the y tick labels
plt.rc('legend', fontsize=16)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
# plt.style.use('tableau-colorblind10')

file_loc = "../data/strain_tests/"
# file_name = "20241031_15_18_25_engagement_time_vs_Fpreload_plotdata"
# file_name = "20241101_15_19_02_engagement_time_vs_Fpreload_plotdata_200hzlowpass_min1to3or50ms"
# file_name = "20241108_14_58_15_timing_vs_Fpreload_plotdata"
file_name = "20241206_23_07_59_timing_vs_Fpreload_plotdata"

data = np.load(save_folder + file_name + ".npy", allow_pickle=True)
# all_V, all_x, all_y, all_yerr, times_by_bucket, forces_by_bucket = data
# all_V, all_x1, all_y1, all_yerr1, all_x2, all_y2, all_yerr2, all_load_cell_engage_times, all_load_cell_disengage_times, all_load_cell_disengage_forces, all_load_cell_preload_forces = data
all_V, all_x1, all_y1, all_yerr1, all_x2, all_y2, all_yerr2, all_load_cell_engage_times, all_load_cell_disengage_times, all_load_cell_disengage_forces, all_load_cell_preload_forces, all_load_cell_initial_forces, all_load_cell_max_forces, all_load_cell_disengage_times_10pct = data

all_freq = sorted(all_load_cell_preload_forces.keys())
num_readings = sum([len(all_load_cell_engage_times[V]) for V in all_V])
num_zeros = sum([len(all_load_cell_engage_times[V][np.where(all_load_cell_engage_times[V] == 0)]) for V in all_V])
num_over_30 = sum([len(all_load_cell_engage_times[V][np.where(all_load_cell_engage_times[V] > 30)]) for V in all_V])
for V in all_V:
    for val in all_load_cell_engage_times[V]:
        if val > 30:
            print("V = {}, val = {}".format(V, val))
print("Num Readings:", num_readings, "Num Zeros:", num_zeros, "Num Over 30:", num_over_30)

# Plot data
colors = ["#1964B0", "#F1932D", "#4DB264", "#DB060B"]  # "#F7F057" = yellow
# fig, axs_all = plt.subplots(2, 2, layout='constrained', figsize=(10, 8), sharex='all', sharey='all')
fig, axs_all = plt.subplots(2, 2, layout='constrained', sharex='all', sharey='all',
                            figsize=(6, 4))
axs = [axs_all[0][0], axs_all[0][1], axs_all[1][0], axs_all[1][1]]
print("Average force for bucket 0", np.mean([xi[0] for xi in all_x1]))

min_x, max_x = np.infty, 0
for i, V in enumerate(all_V):
    min_x = min(min_x, all_x1[i][0])
    max_x = max(max_x, all_x1[i][-1])

for i, V in enumerate(all_V):
    x = all_x1[i]
    y = all_y1[i]
    yerr = all_yerr1[i]
    errorbar = axs[i].errorbar(x, y, yerr=yerr,  # xerr=[np.std(fi[i]) for fi in forces_by_bucket],
                               capsize=5, ecolor='k', elinewidth=2, capthick=2, c='k')  # , alpha=0.4
    scatter = axs[i].scatter(x, y, c='k', s=50)
    # x = all_load_cell_preload_forces[V]
    # y = all_load_cell_engage_times[V]
    # yerr = []
    # axs[i].scatter(x, y, marker='x', c='k', s=50)

    Fpreload_range = np.append(np.linspace(0, max(x) * 0.8, 400), np.linspace(max(x) * 0.8, max_x, 2))
    Fpreload_sim = []
    tengage_sim = []
    for Fpreload in Fpreload_range:
        L_dielectric = 55.5e-3
        w_pin = 2e-3
        depth_pin = 2e-3
        k_dielectric = 54.2
        t_dielectric = 24e-6
        frequency = 1000
        T, Y, dY, X, dX, events, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic = sim_engagement_time(L_dielectric=L_dielectric,
                                                                                                               w_pin=w_pin, depth_pin=depth_pin,
                                                                                                               V_max=V, t_dielectric=t_dielectric,
                                                                                                               period=1 / 2 / frequency, Fext=Fpreload,
                                                                                                               k_dielectric=k_dielectric)
        if type(events) == np.float64:
            Fpreload_sim.append(Fpreload)
            tengage_sim.append(events)
    axs[i].plot(Fpreload_sim, tengage_sim, ls='--', c='tab:red', lw=2)

    # bbox = dict(boxstyle='round', fc='blanchedalmond', ec='orange', alpha=0.5)
    bbox = dict(alpha=0.8, boxstyle='round', fc='white', ec='0.8')
    # axs[i].text(0.95, 0.07, "{} V".format(int(V)), bbox=bbox, transform=axs[i].transAxes, horizontalalignment='right')
    if i == 0 or i == 1 or i == 3:
        axs[i].text(0.05, 0.925, "{} V".format(int(V)), bbox=bbox, transform=axs[i].transAxes, horizontalalignment='left', verticalalignment='top')
    else:
        axs[i].text(0.95, 0.925, "{} V".format(int(V)), bbox=bbox, transform=axs[i].transAxes, horizontalalignment='right', verticalalignment='top')
    axs[i].grid(True)
    # axs[i].set_xlabel("Preload Normal Force (N)")
    # axs[i].set_ylabel("Engagement Time (ms)")
    axs[i].tick_params(labelleft=True)
    # axs[i].set_ylabel("Engagement Time (µs)")
    axs[i].set_xticks([0, 0.5, 1, 1.5])
    # axs[i].set_yticks([-20, 0, 20, 40])
    # axs[i].set_yticks([0, 10, 20])
    axs[i].xaxis.set_tick_params(labelbottom=True)

fig.text(0.00, .98, "(d)", transform=fig.transFigure, horizontalalignment='left', verticalalignment='top')
fig.supxlabel("Preload Normal Force (N)", fontsize=16)
fig.supylabel("Engagement Time (µs)", fontsize=16)
plt.savefig(save_folder + timestamp + ".png", dpi=300)
# # plt.savefig("figures/" + timestamp + ".svg")
plt.savefig(save_folder + timestamp + ".pdf")

plt.show()
