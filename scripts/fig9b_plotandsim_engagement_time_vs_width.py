import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit, root_scalar
from scipy.integrate import solve_ivp
from scipy.special import kv
from scipy.stats import linregress
from datetime import datetime
import csv
# from sim_engagement_time_20241121 import *
from sim_engagement_time import *

now = datetime.now()
name_clarifier = "_engagement_time_vs_width"
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
file_name = "20241208_13_49_12_timing_vs_width_plotdata"

data = np.load(file_loc + file_name + ".npy", allow_pickle=True)
all_V, x1, y1, yerr1, x2, y2, yerr2, all_load_cell_engage_times, all_load_cell_preload_forces, all_load_cell_initial_forces, all_load_cell_max_forces, all_load_cell_disengage_times, all_load_cell_disengage_forces, all_load_cell_disengage_times_10pct = data
print("Min Release Time", np.min(y2), "ms")
print("Average Preload Force", [np.mean(all_load_cell_preload_forces[V]) for V in all_V], "N")

all_freq = sorted(all_load_cell_preload_forces.keys())
num_readings = sum([len(all_load_cell_engage_times[freq]) for freq in all_freq])
num_zeros = sum([len(all_load_cell_engage_times[freq][np.where(all_load_cell_engage_times[freq] == 0)]) for freq in all_freq])
num_over_30 = sum([len(all_load_cell_engage_times[freq][np.where(all_load_cell_engage_times[freq] > 30)]) for freq in all_freq])
for freq in all_freq:
    for val in all_load_cell_engage_times[freq]:
        if val > 30:
            print("freq = {} Hz, val = {}".format(freq, val))
print("Num Readings:", num_readings, "Num Zeros:", num_zeros, "Num Over 30:", num_over_30)

# Plot data
colors = ["#1964B0", "#F1932D", "#4DB264", "#DB060B"]  # "#F7F057" = yellow
fig2, ax2 = plt.subplots(1, 1, layout='constrained', figsize=(6, 3.25))
errorbar2 = ax2.errorbar(x1, y1, yerr=yerr1, capsize=5, ecolor='k', elinewidth=2, capthick=2, color='k', lw=2)
slope, intercept, r, p, se = linregress(x1, y1)
print("Linear fit average data (what's plotted):", slope, intercept, r * r)
slope, intercept, r, p, se = linregress(np.hstack([[width] * len(all_load_cell_engage_times[width]) for width in all_freq]),
                                        np.hstack([all_load_cell_engage_times[width] for width in all_freq]))
print("Linear fit raw data:", slope, intercept, r * r)

scatter2 = ax2.scatter(x1, y1, c='k', s=50)
ax2.grid(True)
ax2.set_xlabel(r"Substrate Width $w_s$ (mm)")
ax2.set_ylabel("Engagement Time (µs)", y=0.4)
bbox = dict(alpha=0.8, boxstyle='round', fc='white', ec='0.8')
# ax2.text(0.96, 0.95, r"$w_s$" + " = 2 mm\nf = 1000 Hz", bbox=bbox, transform=ax2.transAxes, horizontalalignment='right', verticalalignment='top')

w_pin_range = np.linspace(2e-3, 6e-3, 50)
w_sim = []
tengage_sim = []
all_preload_forces = np.hstack([all_load_cell_preload_forces[V] for V in all_V])
Fpreload = np.mean(all_preload_forces)
print("Average preload force:", Fpreload)
for w_pin in w_pin_range:
    L_dielectric = 55.5e-3
    V = 300
    depth_pin = 2e-3
    k_dielectric = 54.2
    t_dielectric = 24e-6
    frequency = 1000
    period = 1 / 2 / frequency
    Fpreload = np.interp(w_pin, all_V, [np.mean(all_load_cell_preload_forces[V]) for V in all_V])  # np.mean(all_load_cell_preload_forces[w_pin * 1e3])
    T, Y, dY, X, dX, events, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic = sim_engagement_time(L_dielectric=L_dielectric,
                                                                                                           w_pin=w_pin, depth_pin=depth_pin,
                                                                                                           V_max=V, t_dielectric=t_dielectric,
                                                                                                           k_dielectric=54.2,
                                                                                                           period=period, Fext=Fpreload)
    if type(events) == np.float64:
        w_sim.append(w_pin * 1e3)
        tengage_sim.append(events)
line_model, = ax2.plot(w_sim, tengage_sim, ls='--', c=colors[-1], lw=2)
ax2.legend([(scatter2, errorbar2), line_model], ["Data", "Model"])
slope, intercept, r, p, se = linregress(w_sim, tengage_sim)
print("Linear fit model:", slope, intercept, r * r)
ax2.text(0.01, .955, "(b)", transform=fig2.transFigure, horizontalalignment='left', verticalalignment='top')

# plt.savefig(save_folder + timestamp + ".png", dpi=300)
# # # plt.savefig("figures/" + timestamp + ".svg")
# plt.savefig(save_folder + timestamp + ".pdf")

plt.show()
