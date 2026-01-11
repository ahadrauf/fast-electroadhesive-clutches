import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

now = datetime.now()
name_clarifier = "_release_time_vs_voltage_10pct"
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

file_name = "20241214_05_46_27_timing_vs_voltage_tr10pct_plotdata"
data = np.load("../data/" + file_name + ".npy", allow_pickle=True)
all_V, x1, y1, yerr1, x2, y2, yerr2, all_load_cell_engage_times, all_load_cell_preload_forces, all_load_cell_initial_forces, all_load_cell_max_forces, all_load_cell_disengage_times, all_load_cell_disengage_forces, all_load_cell_disengage_times_10pct = data
x3 = all_V
y3 = [np.mean(all_load_cell_disengage_times_10pct[width]) for width in all_V]
yerr3 = [np.std(all_load_cell_disengage_times_10pct[width]) for width in all_V]
fig3, ax3 = plt.subplots(1, 1, layout='constrained', figsize=(6, 4))
errorbar3 = ax3.errorbar(x3, y3, yerr=yerr3, capsize=5, ecolor='k', elinewidth=2, capthick=2, color='k', lw=2, ls='--')
scatter3 = ax3.scatter(x3, y3, c='k', s=50, marker='x')

errorbar4 = ax3.errorbar(x2, y2, yerr=yerr3, capsize=5, ecolor='k', elinewidth=2, capthick=2, color='k', lw=2)
scatter4 = ax3.scatter(x2, y2, c='k', s=50, marker='o')
ax3.grid(True)
ax3.set_xlabel("Voltage (V)")
ax3.set_ylabel(r"Load Cell Fall Times (ms)")
ax3.legend([(scatter4, errorbar4), (scatter3, errorbar3)], ["90% Fall Time", "10% Fall Time"])
ax3.text(0.02, .955, "(a)", transform=fig3.transFigure, horizontalalignment='left', verticalalignment='top')

offsets = np.array(y2) - np.array(y3)
print("Average offset between 90% and 10% release times:", offsets, np.mean(offsets), np.std(offsets))

# fig3.savefig(save_folder + timestamp + ".png", dpi=300)
# # # plt.savefig("figures/" + timestamp + ".svg")
# fig3.savefig(save_folder + timestamp + ".pdf")

plt.show()