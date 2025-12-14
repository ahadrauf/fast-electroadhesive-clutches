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

file_name = "20241214_05_46_15_timing_vs_frequency_plotandsim_tr10pct_plotdata"
data = np.load(file_loc + file_name + ".npy", allow_pickle=True)
all_V, all_x1, all_y1, all_yerr1, all_x2, all_y2, all_yerr2, all_load_cell_engage_times, all_load_cell_disengage_times, all_load_cell_disengage_forces, all_load_cell_preload_forces, all_load_cell_initial_forces, all_load_cell_max_forces, all_load_cell_disengage_times_10pct = data
all_freq = sorted(all_load_cell_preload_forces[150].keys())

fig2, axs2_all = plt.subplots(2, 2, layout='constrained', sharex='all', sharey='all', figsize=(6, 4))
axs2 = axs2_all.flat
bbox = dict(alpha=0.8, boxstyle='round', fc='white', ec='0.8')
for i, V in enumerate(all_V):
    x = all_freq
    y = [np.mean(all_load_cell_disengage_times_10pct[V][freq]) for freq in all_freq]
    yerr = [np.std(all_load_cell_disengage_times_10pct[V][freq]) for freq in all_freq]
    errorbar = axs2[i].errorbar(x, y, yerr=yerr, capsize=5, ecolor='k', elinewidth=2, capthick=2, c='k', ls='--')  # alpha=0.4
    scatter = axs2[i].scatter(x, y, c='k', s=50, marker='x')

    x = all_freq
    y = [np.mean(all_load_cell_disengage_times[V][freq]) for freq in all_freq]
    yerr = [np.std(all_load_cell_disengage_times[V][freq]) for freq in all_freq]
    errorbar4 = axs2[i].errorbar(x, y, yerr=yerr, capsize=5, ecolor='k', elinewidth=2, capthick=2, color='k', lw=2)
    scatter4 = axs2[i].scatter(x, y, c='k', s=50, marker='o')

    axs2[i].text(0.95, 0.93, "{} V".format(int(V)), bbox=bbox, transform=axs2[i].transAxes, horizontalalignment='right', verticalalignment='top')
    axs2[i].grid(True)
    axs2[i].tick_params(labelleft=True)
    axs2[i].set_xscale('log')
    axs2[i].set_xticks([0.1, 1, 10, 100, 1000, 10000], ["DC", r"$10^0$", r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$"])
    axs2[i].xaxis.set_tick_params(labelbottom=True)
fig2.supxlabel("Drive Frequency (Hz)", fontsize=16)
fig2.supylabel("Load Cell Fall Times (ms)", fontsize=16)
fig2.text(0.00, .98, "(c)", transform=fig2.transFigure, horizontalalignment='left', verticalalignment='top')

# fig2.savefig(save_folder + timestamp + ".png", dpi=300)
# # # plt.savefig("figures/" + timestamp + ".svg")
# fig2.savefig(save_folder + timestamp + ".pdf")

plt.show()
