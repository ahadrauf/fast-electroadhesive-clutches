import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from parse_data_filenames_20241023 import *
from calculate_timing_stats_20241104 import *
from scipy.stats import zscore

now = datetime.now()
name_clarifier = "_timing_vs_frequency_plotandsim"
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
all_file_names = get_file_dict("20241027_15_25_21", "20241107_15_29_44", key='V', width_filter=w_pin, verbose=False)

all_V = np.array([150, 200, 250, 300])
all_freq = set()
for V, file_names in all_file_names.items():
    for file_name, speed, V, freq, Fpreload, width in file_names:
        all_freq.add(freq)
all_freq = sorted(all_freq)

all_load_cell_engage_times = {V: {freq: [] for freq in all_freq} for V in all_V}
all_load_cell_disengage_forces = {V: {freq: [] for freq in all_freq} for V in all_V}
all_load_cell_disengage_times = {V: {freq: [] for freq in all_freq} for V in all_V}
all_load_cell_preload_forces = {V: {freq: [] for freq in all_freq} for V in all_V}
all_load_cell_initial_forces = {V: {freq: [] for freq in all_freq} for V in all_V}
all_load_cell_max_forces = {V: {freq: [] for freq in all_freq} for V in all_V}
all_load_cell_disengage_times_10pct = {V: {freq: [] for freq in all_freq} for V in all_V}
for V, file_names in all_file_names.items():
    if V not in all_V:
        continue
    for file_name, speed, V, freq, Fpreload, width in file_names:
        print('----------', file_name, '----------')
        depth_pin = get_pin_depth_from_width(width)
        F_preload = calculate_Fpreload(file_loc + file_name, w_pin=width * 1e-3, depth_pin=depth_pin * 1e-3)
        if F_preload > 0.5:
            print("Skipping, Fpreload =", F_preload)
            continue
        release_time, load_cell_before_drop, load_cell_initial, load_cell_max, release_time_10pct = calculate_release_time(file_loc + file_name)
        if release_time > 10:
            print("Skipping, release time unnaturally large:", release_time)
            continue
        engagement_time = calculate_engagement_time(file_loc + file_name)
        if engagement_time > 2000:
            print("Skipping, engagement time unnaturally large:", engagement_time, 'Release Time', release_time)
            continue
        print('Fpreload = {:0.3f}N, Engagement Time = {:0.3f}us, Release Time = {:0.3f}ms, Load Cell Before Drop = {:0.3f}N, Release Time 10pct = {:0.3f}'.format(F_preload, engagement_time, release_time, load_cell_before_drop, release_time_10pct))
        all_load_cell_engage_times[V][freq].append(engagement_time)
        all_load_cell_disengage_times[V][freq].append(release_time)
        all_load_cell_preload_forces[V][freq].append(F_preload)
        all_load_cell_disengage_forces[V][freq].append(load_cell_before_drop)
        all_load_cell_initial_forces[V][freq].append(load_cell_initial)
        all_load_cell_max_forces[V][freq].append(load_cell_max)
        all_load_cell_disengage_times_10pct[V][freq].append(release_time_10pct)

print("------------------ Aggregate Stats Engagement Timing ------------------")
colors = ["#1964B0", "#F1932D", "#4DB264", "#DB060B"]  # "#F7F057" = yellow
all_preload_forces = []
min_disengagement_time = np.infty
for i, V in enumerate(all_V):
    for freq in all_freq:
        all_load_cell_engage_times[V][freq] = np.array(all_load_cell_engage_times[V][freq])
        all_load_cell_disengage_forces[V][freq] = np.array(all_load_cell_disengage_forces[V][freq])
        all_load_cell_disengage_times[V][freq] = np.array(all_load_cell_disengage_times[V][freq])
        all_load_cell_preload_forces[V][freq] = np.array(all_load_cell_preload_forces[V][freq])
        all_load_cell_initial_forces[V][freq] = np.array(all_load_cell_initial_forces[V][freq])
        all_load_cell_max_forces[V][freq] = np.array(all_load_cell_max_forces[V][freq])
        all_load_cell_disengage_times_10pct[V][freq] = np.array(all_load_cell_disengage_times_10pct[V][freq])
        all_preload_forces.extend(all_load_cell_preload_forces[V][freq])

        te_zscores = zscore(all_load_cell_engage_times[V][freq])
        to_delete = all_load_cell_engage_times[V][freq][np.where((te_zscores > 1.5) & (all_load_cell_engage_times[V][freq] > 50))]
        print("Original engagement times:", all_load_cell_engage_times[V][freq])
        print("Zscores:", te_zscores)
        print("To delete:", to_delete)
        for val_to_delete in to_delete:
            all_load_cell_engage_times[V][freq] = all_load_cell_engage_times[V][freq][all_load_cell_engage_times[V][freq] != val_to_delete]

        print("{} V, {} Hz: N = {}, Engage Times = {:0.3f} (std {:0.3f}), raw = {}".format(V, freq, len(all_load_cell_engage_times[V][freq]),
                                                                                           np.mean(all_load_cell_engage_times[V][freq]),
                                                                                           np.std(all_load_cell_engage_times[V][freq]),
                                                                                           all_load_cell_engage_times[V][freq]))
        print("{} V, {} Hz: N = {}, Disengage Times = {:0.3f} (std {:0.3f}), raw = {}".format(V, freq, len(all_load_cell_disengage_times[V][freq]),
                                                                                              np.mean(all_load_cell_disengage_times[V][freq]),
                                                                                              np.std(all_load_cell_disengage_times[V][freq]),
                                                                                              all_load_cell_disengage_times[V][freq]))
        print("{} V, {} Hz: N = {}, Disengage Forces = {:0.3f} (std {:0.3f}), raw = {}".format(V, freq, len(all_load_cell_disengage_forces[V][freq]),
                                                                                               np.mean(all_load_cell_disengage_forces[V][freq]),
                                                                                               np.std(all_load_cell_disengage_forces[V][freq]),
                                                                                               all_load_cell_disengage_forces[V][freq]))
        print("{} V, {} Hz: N = {}, Preload Forces = {:0.3f} (std {:0.3f}), raw = {}".format(V, freq, len(all_load_cell_preload_forces[V][freq]),
                                                                                             np.mean(all_load_cell_preload_forces[V][freq]),
                                                                                             np.std(all_load_cell_preload_forces[V][freq]),
                                                                                             all_load_cell_preload_forces[V][freq]))
        print("{} V, {} Hz: N = {}, Release Times 10% = {:0.3f} (std {:0.3f}), raw = {}".format(V, freq, len(all_load_cell_disengage_times_10pct[V][freq]),
                                                                                                np.mean(all_load_cell_disengage_times_10pct[V][freq]),
                                                                                                np.std(all_load_cell_disengage_times_10pct[V][freq]),
                                                                                                all_load_cell_disengage_times_10pct[V][freq]))

print("Average preload force across all runs", np.nanmean(all_preload_forces), np.nanstd(all_preload_forces))

bbox = dict(alpha=0.8, boxstyle='round', fc='white', ec='0.8')
fig1, ax1_all = plt.subplots(2, 2, figsize=(10, 8), layout='constrained', sharex='all')
axs1 = [ax1_all[0][0], ax1_all[0][1], ax1_all[1][0], ax1_all[1][1]]
all_x1, all_y1, all_yerr1 = [], [], []
for i, V in enumerate(all_V):
    x1 = all_freq
    y1 = [np.mean(all_load_cell_engage_times[V][freq]) for freq in all_freq]
    yerr1 = [np.std(all_load_cell_engage_times[V][freq]) for freq in all_freq]
    errorbar1 = axs1[i].errorbar(x1, y1, yerr=yerr1, capsize=5, ecolor='k', elinewidth=2, capthick=2, color='k', lw=2)
    scatter1 = axs1[i].scatter(x1, y1, c='k', s=50)
    axs1[i].text(0.96, 0.95, "{} V".format(V), bbox=bbox, transform=axs1[i].transAxes, verticalalignment='top', horizontalalignment='right')
    axs1[i].set_xscale('log')
    axs1[i].grid(True)
    axs1[i].set_xlabel("Drive Frequency (Hz)")
    axs1[i].set_ylabel("Engagement Time (µs)")
    axs1[i].set_xticks([0.1, 1, 10, 100, 1000, 10000], ["DC", r"$10^0$", r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$"])

    all_x1.append(x1)
    all_y1.append(y1)
    all_yerr1.append(yerr1)
    print("Minimum engagement time, V = {}: {} us".format(V, np.min(y1)))

fig2, ax2_all = plt.subplots(2, 2, figsize=(10, 8), layout='constrained', sharex='all')
axs2 = [ax2_all[0][0], ax2_all[0][1], ax2_all[1][0], ax2_all[1][1]]
all_x2, all_y2, all_yerr2 = [], [], []
for i, V in enumerate(all_V):
    x2 = all_freq
    y2 = [np.mean(all_load_cell_disengage_times_10pct[V][freq]) for freq in all_freq]
    yerr2 = [np.std(all_load_cell_disengage_times_10pct[V][freq]) for freq in all_freq]
    errorbar2 = axs2[i].errorbar(x2, y2, yerr=yerr2, capsize=5, ecolor='k', elinewidth=2, capthick=2, color='k', lw=2)
    scatter2 = axs2[i].scatter(x2, y2, c='k', s=50)
    axs2[i].text(0.04, 0.95, "{} V".format(V), bbox=bbox, transform=axs2[i].transAxes, verticalalignment='top', horizontalalignment='left')
    axs2[i].set_xscale('log')
    axs2[i].grid(True)
    axs2[i].set_xlabel("Drive Frequency (Hz)")
    axs2[i].set_ylabel("Release Time (ms)")
    axs2[i].set_xticks([0.1, 1, 10, 100, 1000, 10000], ["DC", r"$10^0$", r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$"])

    all_x2.append(x2)
    all_y2.append(y2)
    all_yerr2.append(yerr2)
    print("Minimum release time, V = {}: {} us".format(V, np.min(y2)))
np.save(save_folder + timestamp + "_tr10pct_plotdata", [all_V, all_x1, all_y1, all_yerr1, all_x2, all_y2, all_yerr2, all_load_cell_engage_times,
                                                all_load_cell_disengage_times, all_load_cell_disengage_forces,
                                                all_load_cell_preload_forces, all_load_cell_initial_forces, all_load_cell_max_forces,
                                                all_load_cell_disengage_times_10pct])

# fig2.savefig(save_folder + timestamp + "_release.png", dpi=300)
# # # plt.savefig("figures/" + timestamp + ".svg")
# fig2.savefig(save_folder + timestamp + "_release.pdf")

plt.show()
