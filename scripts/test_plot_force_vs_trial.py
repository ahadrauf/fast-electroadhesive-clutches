from parse_data_filenames import *
from calculate_timing_stats import *
from datetime import datetime

plt.rcParams["font.family"] = "Arial"
plt.rc('font', size=16)  # controls default text size
plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=16)  # fontsize of the y tick labels
plt.rc('legend', fontsize=16)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

now = datetime.now()
name_clarifier = "_plot_force_vs_trial"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
start_time = datetime.now()

save_folder = "../data/"
file_name = "20241206_23_17_24_timing_alldata"
data = np.load("../data/" + file_name + ".npy", allow_pickle=True)
all_V, all_engagement_times, all_release_times, all_release_times_10pct, all_data = np.load(save_folder + file_name + ".npy", allow_pickle=True)

fig, ax = plt.subplots(1, 1, figsize=(6, 4), layout='constrained')
labels = ["DC", r"$10^0$ Hz", r"$10^1$ Hz", r"$10^2$ Hz", r"$10^3$ Hz", r"$10^4$ Hz", r"$10^4$ Hz (Model)"]
ax1_lines = []
markers = ['o', '^', 'd', 'x']
colors = ["#1964B0", "#F1932D", "#4DB264", "#DB060B"]

# file_name, speed, V, freq, Fpreload, width, engagement_time, release_time, load_cell_before_drop, load_cell_initial, load_cell_max, release_time_10pct
dates = np.array([d[0][:17] for d in all_data])
idx = np.argsort(dates)
load_cell_dfs = np.array([d[10] - d[9] for d in all_data])[idx]
# load_cell_dfs_normalized = []
load_cell_dfs_normalized = {}
load_cell_dfs_params = {}
load_cell_dfs_param_counts = {}
for i in range(len(load_cell_dfs)):
    file_name, speed, V, freq, Fpreload, width, engagement_time, release_time, load_cell_before_drop, load_cell_initial, load_cell_max, release_time_10pct = all_data[i]
    # if V != 300 or freq != 1000 or width != 2:
    #     continue
    # if (V, freq, width) != (300.0, 1000.0, 6.0):  # (300, 1000, 2), (250.0, 1000.0, 2.0), (200.0, 1000.0, 2.0), (150.0, 1000.0, 2.0), (300.0, 1000.0, 6.0)
    #     continue
    print(file_name, speed, V, freq, Fpreload, width, engagement_time, release_time, load_cell_before_drop, load_cell_initial, load_cell_max, release_time_10pct)
    if (V, freq, width) not in load_cell_dfs_params.keys():
        load_cell_dfs_params[(V, freq, width)] = load_cell_max - load_cell_initial
        # load_cell_dfs_params[(V, freq, width)] = release_time  # load_cell_max - load_cell_initial
        load_cell_dfs_param_counts[(V, freq, width)] = 0
        load_cell_dfs_normalized[(V, freq, width)] = []
        print("Adding new condition:", V, freq, width)
    # load_cell_dfs_normalized.append((load_cell_max - load_cell_initial) / load_cell_dfs_params[(V, freq, width)])
    load_cell_dfs_normalized[(V, freq, width)].append((load_cell_max - load_cell_initial)/load_cell_dfs_params[(V, freq, width)])
    # load_cell_dfs_normalized.append(release_time)  #  / load_cell_dfs_params[(V, freq, width)])
    load_cell_dfs_param_counts[(V, freq, width)] += 1
print("{} files".format(len(dates)))
for k, v in load_cell_dfs_param_counts.items():
    print(k, "->", v, "times")

min_num_trials = 14

for k in list(load_cell_dfs_normalized.keys()):
    v = load_cell_dfs_normalized[k]
    if len(v) < min_num_trials:
        load_cell_dfs_normalized.pop(k)

for k, v in load_cell_dfs_normalized.items():
    print(k, "After Filtering ->", len(v), "times")

# ax.plot(np.arange(len(load_cell_dfs)), load_cell_dfs)
# ax.plot(np.arange(len(load_cell_dfs_normalized)), load_cell_dfs_normalized)
m, std = np.mean([v[:min_num_trials] for v in load_cell_dfs_normalized.values()], axis=0), np.std([v[:min_num_trials] for v in load_cell_dfs_normalized.values()], axis=0)
ax.plot(np.arange(min_num_trials) + 1, m)
ax.fill_between(np.arange(min_num_trials) + 1, m - std, m + std, color='tab:blue', alpha=0.3)

ax.set_xlabel("Trial Number")
ax.set_ylabel("Shear Force Capacity\nNormalized to First Trial")
# ax.set_xticks([1, 2, 4, 6, 8, 10, 12, 14])
ax.set_xticks(np.arange(1, 15))
ax.axhline(1, color='k', ls='--')
ax.grid(True)

# plt.savefig('../figures/' + timestamp + ".png", dpi=300)
# # # plt.savefig("figures/" + timestamp + ".svg")
# plt.savefig('../figures/' + timestamp + ".pdf")

plt.show()
