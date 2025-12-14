import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit, root_scalar
from scipy.integrate import solve_ivp
from scipy.special import kv
from scipy.stats import linregress
from sim_release_time import *

now = datetime.now()
name_clarifier = "_release_time_vs_voltage"
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
file_name = "20241208_14_27_44_timing_vs_voltage_plotdata"
data = np.load(file_loc + file_name + ".npy", allow_pickle=True)
all_V, x1, y1, yerr1, x2, y2, yerr2, all_load_cell_engage_times, all_load_cell_preload_forces, all_load_cell_initial_forces, all_load_cell_max_forces, all_load_cell_disengage_times, all_load_cell_disengage_forces, all_load_cell_disengage_times_10pct = data

print("Min Release Time", np.min(y2), "ms")
print("Average Preload Force", [np.mean(all_load_cell_preload_forces[V]) for V in all_V], "N")

# Plot data
colors = ["#1964B0", "#F1932D", "#4DB264", "#DB060B"]  # "#F7F057" = yellow
fig2, ax2 = plt.subplots(1, 1, layout='constrained', figsize=(6, 3.25))
# fig2, ax2 = plt.subplots(1, 1, layout='constrained', figsize=(3.45, 1.85))

x2 = all_V
y2 = [np.mean(all_load_cell_disengage_times[width]) for width in all_V]
yerr2 = [np.std(all_load_cell_disengage_times[width]) for width in all_V]

print("Max Time:", np.max(y2))
print("Min Time:", np.min(y2), np.max(y2) / np.min(y2))
print("All Times:", y2)
print("All Err:", yerr2, np.mean(yerr2))

slope, intercept, r, p, se = linregress(x2, y2)
print("Linear fit data:", slope, intercept, r*r)
slope, intercept, r, p, se = linregress(np.hstack([[width]*len(all_load_cell_disengage_times[width]) for width in all_V]),
                                        np.hstack([all_load_cell_disengage_times[width] for width in all_V]))
print("Linear fit raw data:", slope, intercept, r*r)
print(x2)
print(y2)
print(yerr2)

errorbar2 = ax2.errorbar(x2, y2, yerr=yerr2, capsize=5, ecolor='k', elinewidth=2, capthick=2, color='k', lw=2)
scatter2 = ax2.scatter(x2, y2, c='k', s=50)
ax2.grid(True)
ax2.set_xlabel("Voltage (V)")
ax2.set_ylabel("Release Time (ms)")
bbox = dict(alpha=0.8, boxstyle='round', fc='white', ec='0.8')

V_range = np.linspace(100, 300, 3)
V_sim = []
trelease_sim = []
t_release_10pct_sim = []
all_preload_forces = np.hstack([all_load_cell_preload_forces[V] for V in all_V])
all_initial_forces = np.hstack([all_load_cell_initial_forces[V] for V in all_V])
all_disengage_forces = np.hstack([all_load_cell_disengage_forces[V] for V in all_V])
all_max_forces = np.hstack([all_load_cell_max_forces[V] for V in all_V])
# print([all_load_cell_preload_forces[V] for V in all_V])
# print(all_preload_forces)
Fpreload = np.mean(all_preload_forces)
# loadcell_start_ratio = np.mean(np.divide(all_disengage_forces - all_initial_forces, all_max_forces - all_initial_forces))
loadcell_start_ratio = np.mean(np.divide(all_disengage_forces, all_max_forces))
# loadcell_start_ratios_by_voltage = [np.mean(np.divide(all_load_cell_disengage_forces[V] - all_load_cell_initial_forces[V], all_load_cell_max_forces[V] - all_load_cell_initial_forces[V])) for V in all_V]
print("Average preload force:", Fpreload, "std", np.std(all_preload_forces))
print("Average load cell start ratio:", loadcell_start_ratio)
# load_cell_start_ratio = np.mean(all_initial_forces)
if True:
    for V in V_range:
        L_dielectric = 55.5e-3
        w_pin = 2e-3
        depth_pin = 2e-3
        k_dielectric = 54.2
        t_dielectric = 24e-6
        frequency = 1000
        # Fpreload = np.mean(all_load_cell_preload_forces[V])
        # Fpreload = np.interp(V, all_V, [np.mean(all_load_cell_preload_forces[V]) for V in all_V])
        # loadcell_start_ratio = np.interp(V, all_V, loadcell_start_ratios_by_voltage)
        print("Fpreload for V =", V, "=", Fpreload, "start ratio", loadcell_start_ratio)
        T, Y, dY, X, dX, events, t_release_10pct, _, _, _, _, _, _ = sim_release_time(L_dielectric=L_dielectric,
                                                                                      w_pin=w_pin, depth_pin=depth_pin,
                                                                                      V_max=V, t_dielectric=t_dielectric,
                                                                                      frequency=frequency, Fext=Fpreload,
                                                                                      loadcell_start_ratio=loadcell_start_ratio)
        # T, Y, dY, X, dX, events, t_release_10pct, _, _, _, _, _, _ = sim_release_time(L_dielectric=L_dielectric,
        #                                                                                   w_pin=w_pin, depth_pin=depth_pin,
        #                                                                                   V_max=V, t_dielectric=t_dielectric,
        #                                                                                   frequency=frequency, Fext=Fpreload,
        #                                                                                   loadcell_start_ratio=loadcell_start_ratio)
        if type(events) == np.float64:
            V_sim.append(V)
            trelease_sim.append(events)
            t_release_10pct_sim.append(t_release_10pct)
        else:
            print("Output not a float for some reason:", V, events, t_release_10pct)

    print(V_sim)
    print(trelease_sim)
    line_model, = ax2.plot(V_sim, trelease_sim, ls='--', c=colors[-1], lw=2)
    ax2.legend([(scatter2, errorbar2), line_model], ["Data", "Model"])
    slope, intercept, r, p, se = linregress(V_sim, trelease_sim)
    print("Linear fit model:", slope, intercept, r*r)
ax2.text(0.007, .955, "(a)", transform=fig2.transFigure, horizontalalignment='left', verticalalignment='top')

# file_name = "20241214_05_46_27_timing_vs_voltage_tr10pct_plotdata"
# data = np.load("../data/" + file_name + ".npy", allow_pickle=True)
# all_V, x1, y1, yerr1, x2, y2, yerr2, all_load_cell_engage_times, all_load_cell_preload_forces, all_load_cell_initial_forces, all_load_cell_max_forces, all_load_cell_disengage_times, all_load_cell_disengage_forces, all_load_cell_disengage_times_10pct = data
# x3 = all_V
# y3 = [np.mean(all_load_cell_disengage_times_10pct[width]) for width in all_V]
# yerr3 = [np.std(all_load_cell_disengage_times_10pct[width]) for width in all_V]
# fig3, ax3 = plt.subplots(1, 1, layout='constrained', figsize=(6, 4))
# errorbar3 = ax3.errorbar(x3, y3, yerr=yerr3, capsize=5, ecolor='k', elinewidth=2, capthick=2, color='k', lw=2, ls='--')
# scatter3 = ax3.scatter(x3, y3, c='k', s=50, marker='x')
#
# errorbar4 = ax3.errorbar(x2, y2, yerr=yerr3, capsize=5, ecolor='k', elinewidth=2, capthick=2, color='k', lw=2)
# scatter4 = ax3.scatter(x2, y2, c='k', s=50, marker='o')
# ax3.grid(True)
# ax3.set_xlabel("Voltage (V)")
# ax3.set_ylabel(r"Load Cell Fall Times (ms)")
# ax3.legend([(scatter4, errorbar4), (scatter3, errorbar3)], ["90% Fall Time", "10% Fall Time"])
# ax3.text(0.0, .955, "(a)", transform=fig3.transFigure, horizontalalignment='left', verticalalignment='top')
#
# offsets = np.array(y2) - np.array(y3)
# print("Average offset between 90% and 10% release times:", offsets, np.mean(offsets), np.std(offsets))

# fig3.savefig(save_folder + timestamp + ".png", dpi=300)
# # # plt.savefig("figures/" + timestamp + ".svg")
# fig3.savefig(save_folder + timestamp + ".pdf")

plt.show()
