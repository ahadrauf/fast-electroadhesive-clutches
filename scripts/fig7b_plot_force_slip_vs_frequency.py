import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from parse_data_filenames import get_file_dict, get_pin_depth_from_width
from calculate_timing_stats import *
from datetime import datetime

now = datetime.now()
name_clarifier = "_slip_vs_frequency"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
save_folder = "../figures/"

plt.rcParams["font.family"] = "Arial"
plt.rc('font', size=20)  # controls default text size
plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=16)  # fontsize of the y tick labels
plt.rc('legend', fontsize=16)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

file_loc = "../data/"

w_pin = 2
length_pattern = 55.5  # mm
all_file_names = get_file_dict("20241027_15_25_21", "20241107_15_29_44", key='V', width_filter=w_pin, verbose=False)
all_V = np.array([150, 200, 250, 300])
all_freq = set()
for V, file_names in all_file_names.items():
    for file_name, speed, V, freq, Fpreload, width in file_names:
        all_freq.add(freq)
all_freq = sorted(all_freq)

file_name = "20241209_21_29_19_slip_vs_frequency_plotdata"

if file_name == "":
    all_load_cell_max_force = {f: {V: [] for V in all_V} for f in all_freq}
    all_load_cell_slip_recovery_times = {f: {V: [] for V in all_V} for f in all_freq}
    all_load_cell_slip_force_df = {f: {V: [] for V in all_V} for f in all_freq}
    all_sample_delays = []
    for V, file_names in all_file_names.items():
        if V not in all_V:
            continue
        for file_name, speed, V, freq, Fpreload, w_pin in file_names:
            depth_pin = get_pin_depth_from_width(w_pin)
            print('----------', file_name, '----------')
            release_time, load_cell_before_drop, load_cell_initial, load_cell_max, release_time_10pct = calculate_release_time(file_loc + file_name)
            if release_time > 10:
                print("Skipping, release time unnaturally large:", release_time)
                continue

            time_loadcell, load_cell, stage, relay_output, idx_stop_moving, idx_relay_input_rise, idx_relay_input_fall = import_data(file_loc + file_name, zero_load_cell=True, filter_freq=None,
                                                                                                                                     convert_to_kPa=False, convert_time_to_seconds=False)

            # sos = butter(5, (48, 52), btype='bandstop', fs=12800, output='sos')
            # load_cell = sosfiltfilt(sos, load_cell)
            sos = butter(2, 250, btype='lowpass', fs=12800, output='sos')
            load_cell = sosfiltfilt(sos, load_cell)

            max_force = np.max(load_cell[idx_relay_input_rise - 100:idx_relay_input_fall + 100])
            all_load_cell_max_force[freq][V].append(max_force)

            slip_dfs = []
            slip_dts = []
            i = idx_relay_input_rise
            while i < idx_stop_moving - 50:
                didx = 50
                if load_cell[i] == np.max(load_cell[i:i + didx]) and load_cell[i] > max_force / 3:
                    curr_max, curr_max_time = load_cell[i], time_loadcell[i]
                    curr_min, curr_min_time = load_cell[i], time_loadcell[i]
                    while np.min(load_cell[i:i + 10]) < curr_min and np.max(load_cell[i:i + 10]) <= curr_max:
                        if i + 10 >= idx_stop_moving:
                            break
                        curr_min = np.min(load_cell[i:i + 10])
                        curr_min_time = time_loadcell[np.argmin(load_cell[i:i + 10])] + time_loadcell[i]
                        i = i + 10
                    slip_dfs.append((curr_max - curr_min) / curr_max * 100)
                    slip_dts.append((curr_min_time - curr_max_time))
                else:
                    i += 1
            all_load_cell_slip_force_df[freq][V].append(np.nanmax(slip_dfs))
            all_load_cell_slip_recovery_times[freq][V].append(np.median(slip_dts))

            print("Slip %: f = {}, V = {}: {} times: {}".format(freq, V, len(slip_dfs), sorted(slip_dfs, reverse=True)))  # all_load_cell_slip_force_df[f][V][-1] / all_load_cell_max_force[f][V][-1] * 100))
            print("Slip Times: f = {}, V = {}: {} times: {}".format(freq, V, len(slip_dts), sorted(slip_dts)))
else:
    data = np.load("../data/" + file_name + ".npy", allow_pickle=True)
    all_V, all_freq, all_load_cell_slip_force_df, all_load_cell_slip_recovery_times, all_load_cell_max_force = data


print("------------------ Aggregate Stats Force ------------------")
print(all_load_cell_max_force)

fig2, ax2 = plt.subplots(1, 1, layout='constrained', figsize=(6, 4))

markers = ['o', '^', 'd', 'x']
colors = ["#1964B0", "#F1932D", "#4DB264", "#DB060B"]
all_lines = []
for i, V in enumerate(sorted(all_V)):
    y = []
    yerr = []
    for f in sorted(all_load_cell_max_force.keys()):
        print("f = {} Hz, {} V: N = {}, Max Slip = {:0.3f} %, Mean = {:0.3f} %, raw = {}".format(f, V, len(all_load_cell_slip_force_df[f][V]), np.max(all_load_cell_slip_force_df[f][V]),
                                                                                                 np.mean(all_load_cell_slip_force_df[f][V]), all_load_cell_slip_force_df[f][V]))
        y.append(np.mean(all_load_cell_slip_force_df[f][V]))
        yerr.append(np.std(all_load_cell_slip_force_df[f][V]))
    x = all_freq
    errorbar = ax2.errorbar(x=x, y=y, yerr=yerr, lw=1, capsize=5, capthick=1, alpha=0.6, color=colors[i])
    scatter = ax2.scatter(x, y, s=50, color=colors[i], marker=markers[i], zorder=99)
    all_lines.append((errorbar, scatter))

ax2.grid(True)
ax2.set_xlabel("Drive Frequency (Hz)")
ax2.set_ylabel("Max % Shear Force Lost After a Slip", y=0.45)
ax2.yaxis.set_major_formatter(mticker.PercentFormatter())
ax2.set_xscale('log')
ax2.set_xticks([0.1, 1, 10, 100, 1000, 10000], ["DC", r"$10^0$", r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$"])
ax2.set_yticks([0, 25, 50, 75])
ax2.legend(reversed(all_lines), ["300 V", "250 V", "200 V", "150 V"], loc='upper right', fontsize=14)
fig2.text(0.00, 1, "(b)", transform=fig2.transFigure, horizontalalignment='left', verticalalignment='top')

# np.save(save_folder + timestamp + "_plotdata", [all_V, all_freq, all_load_cell_slip_force_df, all_load_cell_slip_recovery_times, all_load_cell_max_force])
# plt.savefig(save_folder + timestamp + ".png", dpi=300)
# # # plt.savefig("figures/" + timestamp + ".svg")
# plt.savefig(save_folder + timestamp + ".pdf")
plt.show()
