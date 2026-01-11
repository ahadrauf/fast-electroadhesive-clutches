import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit, root_scalar
from scipy.integrate import solve_ivp
from scipy.stats import linregress
from datetime import datetime
from sim_release_time import *

now = datetime.now()
name_clarifier = "_release_time_vs_width"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
save_folder = "../data/"

plt.rcParams["font.family"] = "Arial"
plt.rc('font', size=16)  # controls default text size
plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=16)  # fontsize of the y tick labels
plt.rc('legend', fontsize=16)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

file_loc = "../data/strain_tests/"
file_name = "20241208_14_36_48_timing_vs_width_plotdata"
# file_name = "20241214_05_45_49_timing_vs_width_tr10pct_plotdata"

data = np.load(save_folder + file_name + ".npy", allow_pickle=True)
all_widths, x1, y1, yerr1, x2, y2, yerr2, all_load_cell_engage_times, all_load_cell_disengage_times, all_load_cell_disengage_forces, all_load_cell_preload_forces, all_load_cell_initial_forces, all_load_cell_max_forces, all_load_cell_disengage_times_10pct = data
print("Min Release Time", np.min(y2), "ms")

print("Average Preload Force", [np.mean(all_load_cell_preload_forces[V]) for V in all_widths], "N")

# Plot data
colors = ["#1964B0", "#F1932D", "#4DB264", "#DB060B"]
x2 = all_widths
y2 = [np.mean(all_load_cell_disengage_times[width]) for width in all_widths]
yerr2 = [np.std(all_load_cell_disengage_times[width]) for width in all_widths]

slope, intercept, r, p, se = linregress(x2, y2)
print("Linear fit data:", slope, intercept, r * r)
slope, intercept, r, p, se = linregress(np.hstack([[width] * len(all_load_cell_disengage_times[width]) for width in all_widths]),
                                        np.hstack([all_load_cell_disengage_times[width] for width in all_widths]))
print("Linear fit raw data:", slope, intercept, r * r)
print(x2)
print(y2)

fig2, ax2 = plt.subplots(1, 1, layout='constrained', figsize=(6, 3.25))
errorbar2 = ax2.errorbar(x2, y2, yerr=yerr2, capsize=5, ecolor='k', elinewidth=2, capthick=2, color='k', lw=2)
scatter2 = ax2.scatter(x2, y2, c='k', s=50)
ax2.grid(True)
ax2.set_xlabel(r"Substrate Width $w_s$ (mm)")
ax2.set_ylabel("Release Time (ms)")
bbox = dict(alpha=0.8, boxstyle='round', fc='white', ec='0.8')

w_pin_range = np.linspace(2e-3, 6e-3, 3)
w_sim = []
tengage_sim = []
t_release_10pct_sim = []
all_preload_forces = np.hstack([all_load_cell_preload_forces[V] for V in all_widths])
all_initial_forces = np.hstack([all_load_cell_initial_forces[V] for V in all_widths])
all_disengage_forces = np.hstack([all_load_cell_disengage_forces[V] for V in all_widths])
all_max_forces = np.hstack([all_load_cell_max_forces[V] for V in all_widths])
Fpreload = np.mean(all_preload_forces)
loadcell_start_ratio = np.mean(np.divide(all_disengage_forces - all_initial_forces, all_max_forces - all_initial_forces))

print("Average preload force:", Fpreload, "load cell start ratio", loadcell_start_ratio)
print(all_initial_forces)
print(all_max_forces)
if True:
    for w_pin in w_pin_range:
        L_dielectric = 55.5e-3
        V = 300
        depth_pin = 2e-3
        k_dielectric = 54.2
        t_dielectric = 24e-6
        frequency = 1000
        period = 1 / 2 / frequency
        T, Y, dY, X, dX, events, t_release_10pct, _, _, _, _, _, _ = sim_release_time(L_dielectric=L_dielectric,
                                                                                      w_pin=w_pin, depth_pin=depth_pin,
                                                                                      V_max=V, t_dielectric=t_dielectric,
                                                                                      frequency=frequency, Fext=Fpreload,
                                                                                      loadcell_start_ratio=loadcell_start_ratio)
        if type(events) == np.float64:
            w_sim.append(w_pin * 1e3)
            tengage_sim.append(events)
            t_release_10pct_sim.append(t_release_10pct)
        else:
            print("Output not a float for some reason:", V, events, t_release_10pct)

    print(tengage_sim)
    print(t_release_10pct_sim)
    line_model, = ax2.plot(w_sim, tengage_sim, ls='--', c=colors[-1], lw=2)
    ax2.legend([(scatter2, errorbar2), line_model], ["Data", "Model"])

# file_name = "20241214_05_45_49_timing_vs_width_tr10pct_plotdata"
# data = np.load(save_folder + file_name + ".npy", allow_pickle=True)
# all_load_cell_disengage_times_10pct = data[-1]
# # all_V, x1, y1, yerr1, x2, y2, yerr2, all_load_cell_engage_times, all_load_cell_disengage_times, all_load_cell_disengage_forces, all_load_cell_preload_forces, all_load_cell_initial_forces, all_load_cell_max_forces = data
# # all_widths, x1, y1, yerr1, x2, y2, yerr2, all_load_cell_engage_times, all_load_cell_disengage_times, all_load_cell_disengage_forces, all_load_cell_preload_forces, all_load_cell_initial_forces, all_load_cell_max_forces, all_load_cell_disengage_times_10pct = data
# x3 = all_widths
# y3 = [np.mean(all_load_cell_disengage_times_10pct[width]) for width in all_widths]
# yerr3 = [np.std(all_load_cell_disengage_times_10pct[width]) for width in all_widths]
# fig3, ax3 = plt.subplots(1, 1, layout='constrained', figsize=(6, 4))
# errorbar3 = ax3.errorbar(x3, y3, yerr=yerr3, capsize=5, ecolor='k', elinewidth=2, capthick=2, color='k', lw=2, ls='--')
# scatter3 = ax3.scatter(x3, y3, c='k', s=50, marker='x')
# errorbar4 = ax3.errorbar(x2, y2, yerr=yerr3, capsize=5, ecolor='k', elinewidth=2, capthick=2, color='k', lw=2)
# scatter4 = ax3.scatter(x2, y2, c='k', s=50, marker='o')
# ax3.grid(True)
# ax3.set_xlabel(r"Pin Width $w_s$ (mm)")
# # ax3.set_ylabel(r"Release Time $t_{r,90\%}$ (ms)", y=0.4)
# ax3.set_ylabel(r"Load Cell Fall Times (ms)", y=0.5)
# # line_model, = ax3.plot(V_sim, t_release_10pct_sim, ls='--', c=colors[-1], lw=2)
# # ax3.legend([(scatter3, errorbar3), line_model], ["Data", "Model"])
# ax3.legend([(scatter4, errorbar4), (scatter3, errorbar3)], ["90% Fall Time", "10% Fall Time"])
# ax3.text(0.0, .955, "(b)", transform=fig3.transFigure, horizontalalignment='left', verticalalignment='top')
#
# offsets = np.array(y2) - np.array(y3)
# print("Average offset between 90% and 10% release times:", np.mean(offsets), np.std(offsets))

# x2 = all_widths
# y2 = [np.mean(all_load_cell_disengage_times_10pct[width]) for width in all_widths]
# yerr2 = [np.std(all_load_cell_disengage_times_10pct[width]) for width in all_widths]
# fig3, ax3 = plt.subplots(1, 1, layout='constrained', figsize=(3.45, 1.85))
# errorbar3 = ax3.errorbar(x2, y2, yerr=yerr2, capsize=5, ecolor='k', elinewidth=2, capthick=2, color='k', lw=2)
# scatter3 = ax3.scatter(x2, y2, c='k', s=50)
# ax3.grid(True)
# ax3.set_xlabel(r"Pin Width $w_s$ (mm)")
# ax3.set_ylabel(r"Release Time ($t_{r,90\%}$) (ms)")
# line_model, = ax3.plot(w_sim, t_release_10pct_sim, ls='--', c=colors[-1], lw=2)
# ax3.legend([(scatter3, errorbar3), line_model], ["Data", "Model"])
# ax2.text(0.0, .959, "(b)", transform=fig2.transFigure, horizontalalignment='left', verticalalignment='top')
# slope, intercept, r, p, se = linregress(w_sim, tengage_sim)
# print("Linear fit model:", slope, intercept, r*r)

# plt.savefig(save_folder + timestamp + ".png", dpi=300)
# # # plt.savefig("figures/" + timestamp + ".svg")
# plt.savefig(save_folder + timestamp + ".pdf")
# fig3.savefig(save_folder + timestamp + ".png", dpi=300)
# # # plt.savefig("figures/" + timestamp + ".svg")
# fig3.savefig(save_folder + timestamp + ".pdf")

plt.show()
