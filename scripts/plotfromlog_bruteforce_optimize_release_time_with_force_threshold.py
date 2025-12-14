import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

now = datetime.now()
name_clarifier = "_optimize_release_time_with_force_threshold_vs_V"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
save_folder = "../data/"

plt.rcParams["font.family"] = "Arial"
plt.rc('font', size=16)  # controls default text size
plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=16)  # fontsize of the y tick labels
plt.rc('legend', fontsize=16)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# fig, ax = plt.subplots(1, 1, figsize=(6, 5), layout='constrained')
# ax1_right = ax.twinx()
fig, axs = plt.subplots(1, 2, figsize=(12, 4), layout='constrained')
ax, ax1_right = axs

# filenames = ["20251121_18_00_14_optimize_release_time_with_force_threshold",
#              "20251124_00_57_53_optimize_release_time_with_force_threshold",
#              "20251125_15_13_07_optimize_release_time_with_force_threshold"]
# filenames = ["20251125_15_20_15_optimize_release_time_with_force_threshold_wpin=3.0"]
# all_filenames = {2: ["20251121_18_00_14_optimize_release_time_with_force_threshold",
#                      "20251124_00_57_53_optimize_release_time_with_force_threshold",
#                      "20251125_15_13_07_optimize_release_time_with_force_threshold"],
#                  3: ["20251125_15_20_15_optimize_release_time_with_force_threshold_wpin=3.0",
#                      "20251125_23_21_49_optimize_release_time_with_force_threshold_wpin=3.0"],
#                  4: ["20251128_14_50_36_optimize_release_time_with_force_threshold_wpin=4.0",
#                      "20251128_14_51_05_optimize_release_time_with_force_threshold_wpin=4.0",
#                      "20251129_15_51_46_optimize_release_time_with_force_threshold_wpin=4.0"],
#                  5: ["20251129_20_55_03_optimize_release_time_with_force_threshold_wpin=5.0",
#                      "20251129_15_55_08_optimize_release_time_with_force_threshold_wpin=5.0",
#                      "20251129_15_55_50_optimize_release_time_with_force_threshold_wpin=5.0"]}
# 150: ["20251130_22_09_11_optimize_release_time_with_force_threshold_V=150"]
# all_filenames = {200: ["20251130_12_09_16_optimize_release_time_with_force_threshold_V=200"],
#                  250: ["20251130_12_01_53_optimize_release_time_with_force_threshold_V=250"],
#                  300: ["20251130_12_08_46_optimize_release_time_with_force_threshold_V=300"],
#                  350: ["20251201_10_17_50_optimize_release_time_with_force_threshold_V=350"],
#                  400: ["20251130_22_08_38_optimize_release_time_with_force_threshold_V=400"]}
all_filenames = {30: ["20251203_00_52_44_optimize_release_time_with_force_threshold_Ld=30.0"],
                 40: ["20251203_00_56_54_optimize_release_time_with_force_threshold_Ld=40.0"],
                 50: ["20251203_00_49_59_optimize_release_time_with_force_threshold_Ld=50.0"],
                 60: ["20251203_00_53_06_optimize_release_time_with_force_threshold_Ld=60.0"],
                 70: ["20251203_00_57_49_optimize_release_time_with_force_threshold_Ld=70.0"],
                 80: ["20251203_00_58_16_optimize_release_time_with_force_threshold_Ld=80.0"]}
# colors = ['#004488', '#BB5566', '#DDAA33']
# colors = ['#0077BB', '#009988', '#EE7733', '#CC3311']  # '#EE7733', '#0077BB', '#33BBEE', '#EE3377', '#CC3311', '#009988', '#BBBBBB'
# colors = ["#882E72", "#1965B0", "#4EB265", "#E8601C", '#DC050C']
colors = ["#882E72", "#1965B0", "#7BAFDE", "#4EB265", "#E8601C", '#DC050C']

lines = []
lines_right = []
labels = []
labels_right = []
Fshear_min_all = np.linspace(0, 3, 25)
for idx, V in enumerate(all_filenames.keys()):
    filenames = all_filenames[V]
    xopt_all, t_release_all, Fshear_all = [], [], []
    for idy, filename in enumerate(filenames):
        data_all = np.load(save_folder + filename + ".npy", allow_pickle=True).item()
        xopt_all.extend(data_all['xopt_all'])
        t_release_all.extend(data_all['t_release_all'])
        Fshear_all.extend(data_all['Fshear_all'])
        print(V, filename, "F_shearmin/max", np.min(data_all['Fshear_all']), np.max(data_all['Fshear_all']))

    Fshear_best = []
    t_release_best = []
    xopt_best = []
    for Fshear in Fshear_min_all:
        idy = np.where(Fshear_all > Fshear)[0]
        if len(idy) > 0:
            idy = idy[0]
            Fshear_best.append(Fshear)
            t_release_best.append(t_release_all[idy])
            xopt_best.append(xopt_all[idy])

    line_ax, = ax.plot(Fshear_best, t_release_best, color=colors[idx], lw=2)
    line_right, = ax1_right.plot(Fshear_best, xopt_best, color=colors[idx], lw=2)  # , ls='--')
    # ax1_right.plot(Fshear_min_all, Fshear_all, color='tab:red', lw=2)

    lines.append(line_ax)
    lines_right.append(line_right)
    # labels.append(r"$t_r, w_{pin} = $" + str(w_pin) + " mm")
    # labels_right.append(r"$L_{ea}, w_{pin} = $" + str(w_pin) + " mm")
    labels.append("{} V".format(V))
    labels_right.append("{} V".format(V))

# ax1_right.tick_params(axis='y', colors='tab:red')
# ax1_right.yaxis.label.set_color('tab:red')

ax.set_xlabel("Desired Shear Force Capacity (N)")
ax.set_ylabel("Release Time (ms)")
# ax1_right.set_ylabel("Shear Force (N)")
ax1_right.set_xlabel("Desired Shear Force Capacity (N)")
# ax1_right.set_ylabel("Required EA Clutch Length (mm)")
ax.grid(True)
ax1_right.grid(True)

# ax1_right.set_ylim(0, 500)

# lines = list(reversed(lines))
# labels = list(reversed(labels))
ax.legend(lines, labels, fontsize=14, loc='upper left', title='Drive Voltage', title_fontsize=14)
ax1_right.legend(lines_right, labels_right, fontsize=14, loc='upper left', title='Drive Voltage', title_fontsize=14)

# ax.legend(lines + [lines_right[-1]], labels + [labels_right[-1]], fontsize=14, loc='lower right')

fig.text(0, 1, "(b.i)", ha='left', va='top')
# fig.text(0.5, 1, "(b)", ha='left', va='top')

ax1_right.set_ylabel("`", color='white')
fig.text(0.51, 1, '(b.ii)', fontsize=16, ha='left', va='top')
fig.text(0.51, 0.5, r'Required EA Clutch Length $L_s$ (mm)', fontsize=16, ha='left', va='center', rotation=90)
ax1_right.yaxis.set_label_coords(-0.15, 0.5)

bbox = dict(alpha=0.8, boxstyle='round', fc='white', ec='0.8')
fig.text(0.98, 0.04, r"$w_s$ = 2 mm", bbox=bbox, transform=ax.transAxes, horizontalalignment='right', verticalalignment='bottom')
fig.text(0.98, 0.04, r"$w_s$ = 2 mm", bbox=bbox, transform=ax1_right.transAxes, horizontalalignment='right', verticalalignment='bottom')

# axs[i].text(0.95, 0.07, "{} V".format(int(V)), bbox=bbox, transform=axs[i].transAxes, horizontalalignment='right')
# axs[i].text(0.95, 0.93, "{} V".format(int(V)), bbox=bbox, transform=axs[i].transAxes, horizontalalignment='right', verticalalignment='top')
# axs[i].text(0.05, 0.925, "{} V".format(int(V)), bbox=bbox, transform=axs[i].transAxes, horizontalalignment='left', verticalalignment='top')

# plt.savefig("../figures/" + timestamp + ".png", dpi=300)
# # plt.savefig(save_folder + timestamp + ".svg")
# plt.savefig("../figures/" + timestamp + ".pdf")

plt.show()
