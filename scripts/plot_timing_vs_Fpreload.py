import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.stats import zscore
from parse_data_filenames_20241023 import *
from calculate_timing_stats_20241104 import *

now = datetime.now()
name_clarifier = "_timing_vs_Fpreload"
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

# all_V = set()
# for V, file_names in all_file_names.items():
#     for file_name, speed, V, freq, Fpreload, width in file_names:
#         all_V.add(V)
# all_V = sorted(all_V)
all_V = np.array([150, 200, 250, 300])

all_load_cell_engage_times = {V: [] for V in all_V}
all_load_cell_disengage_forces = {V: [] for V in all_V}
all_load_cell_disengage_times = {V: [] for V in all_V}
all_load_cell_preload_forces = {V: [] for V in all_V}
all_load_cell_initial_forces = {V: [] for V in all_V}
all_load_cell_max_forces = {V: [] for V in all_V}
all_load_cell_disengage_times_10pct = {V: [] for V in all_V}
for V, file_names in all_file_names.items():
    if V not in all_V:
        continue
    for file_name, speed, V, freq, Fpreload, width in file_names:
        print('----------', file_name, '----------')
        depth_pin = get_pin_depth_from_width(width)
        F_preload = calculate_Fpreload(file_loc + file_name, w_pin=width * 1e-3, depth_pin=depth_pin * 1e-3)
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
all_disengagement_times = []
all_disengagement_forces = []
min_disengagement_time = np.infty
force_buckets = [(-1, 0.25), (0.25, 0.5), (0.5, 0.77), (0.77, 1.1), (1.1, 1.5)]
engagment_times_by_bucket = [[[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []]]
release_times_by_bucket = [[[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []]]
release_times_10pct_by_bucket = [[[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []]]
forces_by_bucket = [[[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []], [[], [], [], []]]
for i, V in enumerate(all_V):
    all_load_cell_engage_times[V] = np.array(all_load_cell_engage_times[V])
    all_load_cell_disengage_forces[V] = np.array(all_load_cell_disengage_forces[V])
    all_load_cell_disengage_times[V] = np.array(all_load_cell_disengage_times[V])
    all_load_cell_preload_forces[V] = np.array(all_load_cell_preload_forces[V])
    all_load_cell_initial_forces[V] = np.array(all_load_cell_initial_forces[V])
    all_load_cell_max_forces[V] = np.array(all_load_cell_max_forces[V])
    all_load_cell_disengage_times_10pct[V] = np.array(all_load_cell_disengage_times_10pct[V])

    for curr_te, curr_f in zip(all_load_cell_engage_times[V], all_load_cell_preload_forces[V]):
        for idx, bucket in enumerate(force_buckets):
            if bucket[0] <= curr_f <= bucket[1]:
                engagment_times_by_bucket[idx][i].append(curr_te)
                forces_by_bucket[idx][i].append(curr_f)

    for curr_tr, curr_tr10pct, curr_f in zip(all_load_cell_disengage_times[V], all_load_cell_disengage_times_10pct[V], all_load_cell_preload_forces[V]):
        for idx, bucket in enumerate(force_buckets):
            if bucket[0] <= curr_f <= bucket[1]:
                release_times_by_bucket[idx][i].append(curr_tr)
                release_times_10pct_by_bucket[idx][i].append(curr_tr10pct)
                forces_by_bucket[idx][i].append(curr_f)

    # for idx, bucket in enumerate(force_buckets):
    #     print(engagment_times_by_bucket[idx][i])
    #     te_zscores = zscore(engagment_times_by_bucket[idx][i])
    #     if len(te_zscores[~np.isnan(te_zscores)]) > 0:
    #         print(te_zscores)
    #         engagment_times_by_bucket[idx][i] = np.array(engagment_times_by_bucket[idx][i])
    #         engagment_times_by_bucket[idx][i] = engagment_times_by_bucket[idx][i][np.abs(te_zscores) < 2]
    te_zscores = zscore(all_load_cell_engage_times[V])
    to_delete = all_load_cell_engage_times[V][np.where((te_zscores > 0.15) & (all_load_cell_engage_times[V] > 50))]
    print("Original times:", all_load_cell_engage_times[V])
    print("Zscores:", te_zscores)
    print("To delete:", to_delete)
    for val_to_delete in to_delete:
        all_load_cell_preload_forces[V] = all_load_cell_preload_forces[V][all_load_cell_engage_times[V] != val_to_delete]
        all_load_cell_disengage_forces[V] = all_load_cell_disengage_forces[V][all_load_cell_engage_times[V] != val_to_delete]
        all_load_cell_disengage_times[V] = all_load_cell_disengage_times[V][all_load_cell_engage_times[V] != val_to_delete]
        all_load_cell_initial_forces[V] = all_load_cell_initial_forces[V][all_load_cell_engage_times[V] != val_to_delete]
        all_load_cell_max_forces[V] = all_load_cell_max_forces[V][all_load_cell_engage_times[V] != val_to_delete]
        all_load_cell_disengage_times_10pct[V] = all_load_cell_disengage_times_10pct[V][all_load_cell_engage_times[V] != val_to_delete]
        all_load_cell_engage_times[V] = all_load_cell_engage_times[V][all_load_cell_engage_times[V] != val_to_delete]

    for idx, bucket in enumerate(force_buckets):
        engagment_times_by_bucket[idx][i] = np.array(engagment_times_by_bucket[idx][i])
        print("Bucket", idx, "Before Engagement Times:", engagment_times_by_bucket[idx][i])
        for val_to_delete in to_delete:
            engagment_times_by_bucket[idx][i] = engagment_times_by_bucket[idx][i][engagment_times_by_bucket[idx][i] != val_to_delete]
        print("Bucket", idx, "After Engagement Times:", engagment_times_by_bucket[idx][i])

    print("{} V: N = {}, Engage Times = {:0.3f} (std {:0.3f}), raw = {}".format(V, len(all_load_cell_engage_times[V]),
                                                                                np.mean(all_load_cell_engage_times[V]),
                                                                                np.std(all_load_cell_engage_times[V]),
                                                                                list(all_load_cell_engage_times[V])))
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
print("Minimum engagement times:", [np.min(all_load_cell_engage_times[V]) for V in all_V])
print("Average preload force across all readings:", np.mean([np.mean(all_load_cell_preload_forces[V]) for V in all_V]),
      np.std([np.mean(all_load_cell_preload_forces[V]) for V in all_V]))

bbox = dict(alpha=0.8, boxstyle='round', fc='white', ec='0.8')
fig1, ax1_all = plt.subplots(2, 2, figsize=(10, 8), layout='constrained', sharex='all')
axs1 = [ax1_all[0][0], ax1_all[0][1], ax1_all[1][0], ax1_all[1][1]]
all_x1, all_y1, all_yerr1 = [], [], []
for i, V in enumerate(all_V):
    x1 = [np.mean(fi[i]) for fi in forces_by_bucket]
    y1 = [np.mean(ti[i]) for ti in engagment_times_by_bucket]
    yerr1 = [np.std(ti[i]) for ti in engagment_times_by_bucket]
    errorbar1 = axs1[i].errorbar(x1, y1, yerr=yerr1, capsize=5, ecolor='k', elinewidth=2, capthick=2, color='k', lw=2)
    scatter1 = axs1[i].scatter(x1, y1, c='k', s=50)
    # x1 = all_load_cell_preload_forces[V]
    # y1 = all_load_cell_engage_times[V]
    # yerr1 = []
    # axs1[i].scatter(x1, y1, marker='x', c='k', s=50)
    axs1[i].text(0.95, 0.07, "{} V".format(V), bbox=bbox, transform=axs1[i].transAxes, horizontalalignment='right')
    axs1[i].grid(True)
    axs1[i].set_xlabel("Preload Normal Force (N)")
    axs1[i].set_ylabel("Engagement Time (µs)")
    axs1[i].set_xticks([0, 0.5, 1, 1.5])

    all_x1.append(x1)
    all_y1.append(y1)
    all_yerr1.append(yerr1)
    print("Minimum engagement time, V = {}: {} us".format(V, np.min(y1)))

fig2, ax2_all = plt.subplots(2, 2, figsize=(10, 8), layout='constrained', sharex='all')
axs2 = [ax2_all[0][0], ax2_all[0][1], ax2_all[1][0], ax2_all[1][1]]
all_x2, all_y2, all_yerr2 = [], [], []
for i, V in enumerate(all_V):
    x2 = [np.mean(fi[i]) for fi in forces_by_bucket]
    y2 = [np.mean(ti[i]) for ti in release_times_by_bucket]
    yerr2 = [np.std(ti[i]) for ti in release_times_by_bucket]
    # y2 = [np.mean(ti[i]) for ti in release_times_10pct_by_bucket]
    # yerr2 = [np.std(ti[i]) for ti in release_times_10pct_by_bucket]
    errorbar2 = axs2[i].errorbar(x2, y2, yerr=yerr2, capsize=5, ecolor='k', elinewidth=2, capthick=2, color='k', lw=2)
    scatter2 = axs2[i].scatter(x2, y2, c='k', s=50)
    # x2 = all_load_cell_preload_forces[V]
    # y2 = all_load_cell_disengage_times[V]
    # yerr2 = []
    # axs2[i].scatter(x2, y2, marker='x', c='k', s=50)
    axs2[i].text(0.95, 0.07, "{} V".format(V), bbox=bbox, transform=axs2[i].transAxes, horizontalalignment='right')
    axs2[i].grid(True)
    axs2[i].set_xlabel("Preload Normal Force (N)")
    axs2[i].set_ylabel("Release Time (ms)")
    axs2[i].set_xticks([0, 0.5, 1, 1.5])

    all_x2.append(x2)
    all_y2.append(y2)
    all_yerr2.append(yerr2)
    print("Minimum release time, V = {}: {} us".format(V, np.min(y2)))
np.save(save_folder + timestamp + "_tr10pct_plotdata", [all_V, all_x1, all_y1, all_yerr1, all_x2, all_y2, all_yerr2, all_load_cell_engage_times,
                                                all_load_cell_disengage_times, all_load_cell_disengage_forces,
                                                all_load_cell_preload_forces, all_load_cell_initial_forces, all_load_cell_max_forces, all_load_cell_disengage_times_10pct])

# fig2.savefig(save_folder + timestamp + "_release.png", dpi=300)
# # # plt.savefig("figures/" + timestamp + ".svg")
# fig2.savefig(save_folder + timestamp + "_release.pdf")

plt.show()
