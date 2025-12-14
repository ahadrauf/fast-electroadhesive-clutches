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
name_clarifier = "_engagement_time_vs_width_bars2025"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
save_folder = "../figures/"

plt.rcParams["font.family"] = "Arial"
plt.rc('font', size=26)  # controls default text size
plt.rc('axes', labelsize=26)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=26)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=26)  # fontsize of the y tick labels
plt.rc('legend', fontsize=20)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
# plt.style.use('tableau-colorblind10')

file_loc = "../data/"
file_name = "20241108_15_52_02_timing_vs_width_plotdata"
# file_name = "20241204_00_17_01_timing_vs_width_plotdata"
file_name = "20241206_23_08_07_timing_vs_width_plotdata"
file_name = "20241208_13_49_12_timing_vs_width_plotdata_Fpreloadunder0.25"

data = np.load(file_loc + file_name + ".npy", allow_pickle=True)
# all_V, x1, y1, yerr1, x2, y2, yerr2, all_load_cell_engage_times, all_load_cell_disengage_times, all_load_cell_disengage_forces, all_load_cell_preload_forces = data
# all_V, x1, y1, yerr1, x2, y2, yerr2, all_load_cell_engage_times, all_load_cell_disengage_times, all_load_cell_disengage_forces, all_load_cell_preload_forces, all_load_cell_initial_forces, all_load_cell_max_forces = data
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
# num_readings = sum([sum([len(all_load_cell_engage_times[V][freq]) for freq in all_freq]) for V in all_V])
# num_zeros = sum([sum([len(all_load_cell_engage_times[V][freq][np.where(all_load_cell_engage_times[V][freq] == 0)]) for freq in all_freq]) for V in all_V])
# num_over_30 = sum([sum([len(all_load_cell_engage_times[V][freq][np.where(all_load_cell_engage_times[V][freq] > 30)]) for freq in all_freq]) for V in all_V])
# for V in all_V:
#     for freq in all_freq:
#         for val in all_load_cell_engage_times[V][freq]:
#             if val > 30:
#                 print("V = {}, freq = {} Hz, val = {}".format(V, freq, val))
print("Num Readings:", num_readings, "Num Zeros:", num_zeros, "Num Over 30:", num_over_30)

# Plot data
colors = ["#1964B0", "#F1932D", "#4DB264", "#DB060B"]  # "#F7F057" = yellow
# fig2, ax2 = plt.subplots(1, 1, layout='constrained', figsize=(6, 3.25))
fig2, ax2 = plt.subplots(1, 1, layout='constrained', figsize=(5.5, 4))
errorbar2 = ax2.errorbar(x1, y1, yerr=yerr1, capsize=5, ecolor='k', elinewidth=2, capthick=2, color='k', lw=2)
slope, intercept, r, p, se = linregress(x1, y1)
print("Linear fit average data (what's plotted):", slope, intercept, r * r)
slope, intercept, r, p, se = linregress(np.hstack([[width] * len(all_load_cell_engage_times[width]) for width in all_freq]),
                                        np.hstack([all_load_cell_engage_times[width] for width in all_freq]))
print("Linear fit raw data:", slope, intercept, r * r)
print(x1)
print(y1)
print(list(np.hstack([[width] * len(all_load_cell_engage_times[width]) for width in all_freq])))
print(list(np.hstack([all_load_cell_engage_times[width] for width in all_freq])))
scatter2 = ax2.scatter(x1, y1, c='k', s=50)
ax2.grid(True)
ax2.set_xlabel(r"Substrate Width (mm)")
ax2.set_ylabel("Engagement Time (µs)", y=0.4)
bbox = dict(alpha=0.8, boxstyle='round', fc='white', ec='0.8')
# ax2.text(0.96, 0.95, r"$w_s$" + " = 2 mm\nf = 1000 Hz", bbox=bbox, transform=ax2.transAxes, horizontalalignment='right', verticalalignment='top')

w_pin_range = np.linspace(2e-3, 6e-3, 50)
w_sim = []
tengage_sim = []
all_preload_forces = np.hstack([all_load_cell_preload_forces[V] for V in all_V])
# print([all_load_cell_preload_forces[V] for V in all_V])
# print(all_preload_forces)
Fpreload = np.mean(all_preload_forces)
# Fpreload = np.max([np.mean(all_load_cell_preload_forces[V]) for V in all_V])
print("Average preload force:", Fpreload)
# for w_pin in w_pin_range:
#     L_dielectric = 55.5e-3
#     V = 300
#     depth_pin = 2e-3
#     k_dielectric = 54.2
#     t_dielectric = 24e-6
#     frequency = 1000
#     period = 1 / 2 / frequency
#     Fpreload = np.interp(w_pin, all_V, [np.mean(all_load_cell_preload_forces[V]) for V in all_V])  # np.mean(all_load_cell_preload_forces[w_pin * 1e3])
#     T, Y, dY, X, dX, events, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic = sim_engagement_time(L_dielectric=L_dielectric,
#                                                                                                            w_pin=w_pin, depth_pin=depth_pin,
#                                                                                                            V_max=V, t_dielectric=t_dielectric,
#                                                                                                            k_dielectric=54.2,
#                                                                                                            period=period, Fext=Fpreload)
#     if type(events) == np.float64:
#         w_sim.append(w_pin * 1e3)
#         tengage_sim.append(events)

w_sim = [2.0, 2.0816326530612246, 2.163265306122449, 2.2448979591836733, 2.326530612244898, 2.4081632653061225, 2.489795918367347, 2.5714285714285716, 2.653061224489796, 2.7346938775510203, 2.8163265306122454, 2.8979591836734695, 2.979591836734694, 3.0612244897959187, 3.142857142857143, 3.2244897959183674, 3.306122448979592, 3.3877551020408165, 3.4693877551020407, 3.5510204081632657, 3.63265306122449, 3.7142857142857144, 3.795918367346939, 3.8775510204081636, 3.9591836734693877, 4.040816326530613, 4.122448979591837, 4.204081632653061, 4.285714285714286, 4.36734693877551, 4.448979591836735, 4.530612244897959, 4.612244897959184, 4.6938775510204085, 4.775510204081633, 4.857142857142857, 4.938775510204081, 5.020408163265306, 5.102040816326531, 5.183673469387755, 5.26530612244898, 5.346938775510203, 5.428571428571429, 5.510204081632653, 5.591836734693878, 5.673469387755102, 5.755102040816326, 5.836734693877551, 5.918367346938775, 6.0]
tengage_sim = [2.2350996434847064, 2.3065412237796266, 2.3786575570629105, 2.451116202705798, 2.524084122232585, 2.597423542235959, 2.6712439078482006, 2.7454369472750275, 2.8201036251702747, 2.8951230292044956, 2.970380051304579, 3.046157448374348, 3.1224193997103846, 3.1888073973811633, 3.2312339921720143, 3.3100043723263717, 3.3892459129509835, 3.4686781951201007, 3.548435912384716, 3.6286393009699465, 3.7090611155556052, 3.7897523964010045, 3.870783560800724, 3.9521294756827814, 4.03380306980721, 4.115670061465796, 4.1979716692587825, 4.2806600011323575, 4.363537624824963, 4.446923891840757, 4.530513809776878, 4.614503288212463, 4.698830463302549, 4.78357112902138, 4.868576856623812, 4.953906314934186, 5.039616773311852, 5.125747038024387, 5.212238105288806, 5.299036933624849, 5.3863021103660005, 5.473976528758725, 5.561965335691823, 5.650341273828851, 5.673716059114549, 5.759287879142215, 5.845247886150517, 5.92677230556257, 6.013192774407166, 6.100094365463364]

print(w_sim)
print(tengage_sim)
line_model, = ax2.plot(w_sim, tengage_sim, ls='--', c=colors[-1], lw=2)
ax2.legend([(scatter2, errorbar2), line_model], ["Data", "Model"], loc='upper left')
# ax2.text(0.005, .96, "(b)", transform=fig2.transFigure, horizontalalignment='left', verticalalignment='top')
slope, intercept, r, p, se = linregress(w_sim, tengage_sim)
print("Linear fit model:", slope, intercept, r * r)
print(ax2.get_ylim())
ax2.set_ylim(-4.646277307275232, 22.201907321532964)
xlim = ax2.get_xlim()
ax2.fill_between(ax2.get_xlim(), y1=ax2.get_ylim()[0], y2=15, color='#FFC000', alpha=0.4)
ax2.set_xlim(xlim)

# (-4.646277307275232, 22.201907321532964)
# (-4.5284140558353085, 17.917407017806095)

# plt.savefig(save_folder + timestamp + ".png", dpi=300)
# plt.savefig(save_folder + timestamp + ".svg")
# plt.savefig(save_folder + timestamp + ".pdf")

plt.show()
