import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from parse_data_filenames_20241023 import *
from calculate_timing_stats_20241104 import *
from scipy.stats import zscore

now = datetime.now()
name_clarifier = "_timing_vs_voltage"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
save_folder = "C:/Users/ahadrauf/Desktop/Research/latex/electroadhesive_dynamics_paper/figures_test/"

plt.rcParams["font.family"] = "Arial"
# plt.rcParams["font.family"] = "serif"
# plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rc('font', size=24)  # controls default text size
plt.rc('axes', labelsize=20)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=20)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=20)  # fontsize of the y tick labels
plt.rc('legend', fontsize=20)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
# plt.style.use('tableau-colorblind10')

file_loc = "../data/strain_tests/"

w_pin = 2
depth_pin = 2
length_pattern = 55.5  # mm
frequency = 1000
all_file_names = get_file_dict("20241027_15_25_21", "20241103_18_59_41", key='V', width_filter=w_pin, f_filter=frequency, verbose=False)

all_V = set()
for V, file_names in all_file_names.items():
    for file_name, speed, V, freq, Fpreload, width in file_names:
        all_V.add(V)
all_V = sorted(all_V)

all_load_cell_engage_times = {V: [] for V in all_V}
all_load_cell_disengage_forces = {V: [] for V in all_V}
all_load_cell_disengage_times = {V: [] for V in all_V}
all_load_cell_preload_forces = {V: [] for V in all_V}
all_load_cell_initial_forces = {V: [] for V in all_V}
all_load_cell_max_forces = {V: [] for V in all_V}
all_load_cell_disengage_times_10pct = {V: [] for V in all_V}
all_engagement_times = []
for V, file_names in all_file_names.items():
    for file_name, speed, V, freq, Fpreload, width in file_names:
        print('----------', file_name, '----------')
        depth_pin = get_pin_depth_from_width(width)
        F_preload = calculate_Fpreload(file_loc + file_name, w_pin=width * 1e-3, depth_pin=depth_pin * 1e-3)
        if F_preload > 0.25:
            print("Skipping, Fpreload =", F_preload)
            continue
        engagement_time = calculate_engagement_time(file_loc + file_name)
        release_time, load_cell_before_drop, load_cell_initial, load_cell_max, release_time_10pct = calculate_release_time(file_loc + file_name)
        if engagement_time > 2000:
            print("Skipping, engagement time unnaturally large:", engagement_time, 'Release Time', release_time)
            continue
        if release_time > 10:
            print("Skipping, release time unnaturally large:", release_time, 'Engagement Time', engagement_time)
            continue
        # if F_preload > 0.25:
        #     print("Skipping, Fpreload =", F_preload, "Engagement Time:", engagement_time, "Release Time:", release_time)
        #     continue
        print('Fpreload = {:0.3f}N, Engagement Time = {:0.3f}us, Release Time = {:0.3f}ms, Load Cell Before Drop = {:0.3f}N, Release Time 10pct = {:0.3f}'.format(F_preload, engagement_time, release_time, load_cell_before_drop, release_time_10pct))
        all_load_cell_engage_times[V].append(engagement_time)
        all_load_cell_disengage_times[V].append(release_time)
        all_load_cell_preload_forces[V].append(F_preload)
        all_load_cell_disengage_forces[V].append(load_cell_before_drop)
        all_load_cell_initial_forces[V].append(load_cell_initial)
        all_load_cell_max_forces[V].append(load_cell_max)
        all_load_cell_disengage_times_10pct[V].append(release_time_10pct)

print("------------------ Aggregate Stats Engagement Timing ------------------")
all_freq = np.array(sorted(all_load_cell_engage_times.keys()))
colors = ["#1964B0", "#F1932D", "#4DB264", "#DB060B"]  # "#F7F057" = yellow
# fig, ax1 = plt.subplots(1, 1, layout='constrained')
fig1, ax1 = plt.subplots(1, 1, layout='constrained')
fig2, ax2 = plt.subplots(1, 1, layout='constrained')
all_disengagement_times = []
all_disengagement_forces = []
min_disengagement_time = np.infty
for i, V in enumerate(all_V):
    all_load_cell_engage_times[V] = np.array(all_load_cell_engage_times[V])
    all_load_cell_disengage_forces[V] = np.array(all_load_cell_disengage_forces[V])
    all_load_cell_disengage_times[V] = np.array(all_load_cell_disengage_times[V])
    all_load_cell_preload_forces[V] = np.array(all_load_cell_preload_forces[V])
    all_load_cell_initial_forces[V] = np.array(all_load_cell_initial_forces[V])
    all_load_cell_max_forces[V] = np.array(all_load_cell_max_forces[V])
    all_load_cell_disengage_times_10pct[V] = np.array(all_load_cell_disengage_times_10pct[V])

    te_zscores = zscore(all_load_cell_engage_times[V])
    to_delete = all_load_cell_engage_times[V][np.where((te_zscores > 0.3) & (all_load_cell_engage_times[V] > 50))]
    print("Original engagement times:", all_load_cell_engage_times[V])
    print("Zscores:", te_zscores)
    print("To delete:", to_delete)
    for val_to_delete in to_delete:
        all_load_cell_engage_times[V] = all_load_cell_engage_times[V][all_load_cell_engage_times[V] != val_to_delete]

    print("{} V: N = {}, Engage Times = {:0.3f} (std {:0.3f}), raw = {}".format(V, len(all_load_cell_engage_times[V]),
                                                                                np.mean(all_load_cell_engage_times[V]),
                                                                                np.std(all_load_cell_engage_times[V]),
                                                                                all_load_cell_engage_times[V]))
    print("{} V: N = {}, Disengage Times = {:0.3f} (std {:0.3f}), raw = {}".format(V, len(all_load_cell_disengage_times[V]),
                                                                                   np.mean(all_load_cell_disengage_times[V]),
                                                                                   np.std(all_load_cell_disengage_times[V]),
                                                                                   all_load_cell_disengage_times[V]))
    print("{} V: N = {}, Disengage Forces = {:0.3f} (std {:0.3f}), raw = {}".format(V, len(all_load_cell_disengage_forces[V]),
                                                                                    np.mean(all_load_cell_disengage_forces[V]),
                                                                                    np.std(all_load_cell_disengage_forces[V]),
                                                                                    all_load_cell_disengage_forces[V]))
    print("{} V: N = {}, Preload Forces = {:0.3f} (std {:0.3f}), raw = {}".format(V, len(all_load_cell_preload_forces[V]),
                                                                                  np.mean(all_load_cell_preload_forces[V]),
                                                                                  np.std(all_load_cell_preload_forces[V]),
                                                                                  all_load_cell_preload_forces[V]))
    print("{} V: N = {}, Release Time 10% = {:0.3f} (std {:0.3f}), raw = {}".format(V, len(all_load_cell_disengage_times_10pct[V]),
                                                                                    np.mean(all_load_cell_disengage_times_10pct[V]),
                                                                                    np.std(all_load_cell_disengage_times_10pct[V]),
                                                                                    all_load_cell_disengage_times_10pct[V]))
print("Minimum engagement times:", [np.min(all_load_cell_engage_times[V]) for V in all_V])
print("Average preload force across all readings:", np.mean([np.mean(all_load_cell_preload_forces[V]) for V in all_V]),
      np.std([np.mean(all_load_cell_preload_forces[V]) for V in all_V]))

x1 = all_V
y1 = [np.mean(all_load_cell_engage_times[V]) for V in all_V]
yerr1 = [np.std(all_load_cell_engage_times[V]) for V in all_V]
errorbar1 = ax1.errorbar(x1, y1, yerr=yerr1, capsize=5, ecolor='k', elinewidth=2, capthick=2, color='k', lw=2)
scatter1 = ax1.scatter(x1, y1, c='k', s=50)
bbox = dict(alpha=0.8, boxstyle='round', fc='white', ec='0.8')
# axs[i].text(0.95, 0.07, "{} V".format(V), bbox=bbox, transform=axs[i].transAxes, horizontalalignment='right')
# ax1[i].text(0.95, 0.95, r"$w_s$ = 2 mm\n1000 Hz", bbox=bbox, transform=ax1[i].transAxes, horizontalalignment='right', verticalalignment='top')
ax1.grid(True)
ax1.set_xlabel("Voltage (V)")
ax1.set_ylabel("Engagement Time (µs)")
# ax1[i].set_xticks([0, 0.5, 1, 1.5])

x2 = all_V
y2 = [np.mean(all_load_cell_disengage_times[V]) for V in all_V]
yerr2 = [np.std(all_load_cell_disengage_times[V]) for V in all_V]
# y2 = [np.mean(all_load_cell_disengage_times_10pct[V]) for V in all_V]
# yerr2 = [np.std(all_load_cell_disengage_times_10pct[V]) for V in all_V]
errorbar2 = ax2.errorbar(x2, y2, yerr=yerr2, capsize=5, ecolor='k', elinewidth=2, capthick=2, color='k', lw=2)
scatter2 = ax2.scatter(x2, y2, c='k', s=50)
ax2.grid(True)
ax2.set_xlabel("Voltage (V)")
ax2.set_ylabel("Release Time (ms)")
np.save(save_folder + timestamp + "_tr10pct_plotdata", [all_V, x1, y1, yerr1, x2, y2, yerr2, all_load_cell_engage_times,
                                                all_load_cell_preload_forces, all_load_cell_initial_forces, all_load_cell_max_forces,
                                                all_load_cell_disengage_times, all_load_cell_disengage_forces,
                                                all_load_cell_disengage_times_10pct])

# fig2.savefig(save_folder + timestamp + "_release.png", dpi=300)
# # # plt.savefig("figures/" + timestamp + ".svg")
# fig2.savefig(save_folder + timestamp + "_release.pdf")

plt.show()
