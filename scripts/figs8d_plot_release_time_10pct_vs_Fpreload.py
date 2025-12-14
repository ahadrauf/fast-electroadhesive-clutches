import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

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

file_name = "20241214_05_46_00_timing_vs_Fpreload_tr10pct_plotdata"
data = np.load("../data/" + file_name + ".npy", allow_pickle=True)
all_V, all_x1, all_y1, all_yerr1, all_x2, all_y2, all_yerr2, all_load_cell_engage_times, all_load_cell_disengage_times, all_load_cell_disengage_forces, all_load_cell_preload_forces, all_load_cell_initial_forces, all_load_cell_max_forces, all_load_cell_disengage_times_10pct = data
fig2, axs2_all = plt.subplots(2, 2, layout='constrained', sharex='all', sharey='all', figsize=(6, 4))
axs2 = axs2_all.flat

force_buckets = [(-1, 0.25), (0.25, 0.5), (0.5, 0.77), (0.77, 1.1), (1.1, 1.5)]
engagment_times_by_bucket = [[[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []]]
release_times_by_bucket = [[[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []]]
release_times_10pct_by_bucket = [[[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []]]
forces_by_bucket = [[[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []]]
bbox = dict(alpha=0.8, boxstyle='round', fc='white', ec='0.8')
for i, V in enumerate(all_V):
    all_load_cell_disengage_times_10pct[V] = np.array(all_load_cell_disengage_times_10pct[V])
    for curr_tr, curr_tr10pct, curr_f in zip(all_load_cell_disengage_times[V], all_load_cell_disengage_times_10pct[V], all_load_cell_preload_forces[V]):
        for idx, bucket in enumerate(force_buckets):
            if bucket[0] <= curr_f <= bucket[1]:
                forces_by_bucket[idx][i].append(curr_f)
                release_times_by_bucket[idx][i].append(curr_tr)
                release_times_10pct_by_bucket[idx][i].append(curr_tr10pct)

for i, V in enumerate(all_V):
    x = [np.mean(fi[i]) for fi in forces_by_bucket]
    y = [np.mean(ti[i]) for ti in release_times_10pct_by_bucket]
    yerr = [np.std(ti[i]) for ti in release_times_10pct_by_bucket]
    errorbar = axs2[i].errorbar(x, y, yerr=yerr, capsize=5, ecolor='k', elinewidth=2, capthick=2, c='k', ls='--')  # alpha=0.4
    scatter = axs2[i].scatter(x, y, c='k', s=50, marker='x')

    x2 = [np.mean(fi[i]) for fi in forces_by_bucket]
    y2 = [np.mean(ti[i]) for ti in release_times_by_bucket]
    yerr2 = [np.std(ti[i]) for ti in release_times_by_bucket]
    errorbar4 = axs2[i].errorbar(x2, y2, yerr=yerr2, capsize=5, ecolor='k', elinewidth=2, capthick=2, color='k', lw=2)
    scatter4 = axs2[i].scatter(x2, y2, c='k', s=50, marker='o')

    axs2[i].text(0.05, 0.93, "{} V".format(int(V)), bbox=bbox, transform=axs2[i].transAxes, va='top', horizontalalignment='left')
    axs2[i].grid(True)
    axs2[i].set_xticks([0, 0.5, 1, 1.5])
    axs2[i].set_yticks([0, 1, 2])
    axs2[i].yaxis.set_tick_params(labelleft=True)
    axs2[i].xaxis.set_tick_params(labelbottom=True)

    offsets = np.array(y2) - np.array(y)
    print("Average offset between 90% and 10% release times at V =", V, np.mean(offsets), np.std(offsets))
fig2.supxlabel("Preload Normal Force (N)", fontsize=16)
fig2.supylabel("Load Cell Fall Times (ms)", fontsize=16)
fig2.text(0.00, .98, "(d)", transform=fig2.transFigure, horizontalalignment='left', verticalalignment='top')

# fig2.savefig(save_folder + timestamp + ".png", dpi=300)
# # # plt.savefig("figures/" + timestamp + ".svg")
# fig2.savefig(save_folder + timestamp + ".pdf")

plt.show()
