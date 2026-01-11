import numpy as np
import matplotlib.pyplot as plt
from import_strain_test_data import import_data
from datetime import datetime
from scipy.signal import savgol_filter, butter, sosfiltfilt
from scipy.stats import linregress

now = datetime.now()
name_clarifier = "_force_vs_time_sample"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
save_folder = "C:/Users/ahadrauf/Desktop/Research/latex/electroadhesive_dynamics_paper/figures_test/"

plt.rcParams["font.family"] = "Arial"
plt.rc('font', size=16)  # controls default text size
plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=16)  # fontsize of the y tick labels
plt.rc('legend', fontsize=16)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

file_loc = "C:/Users/ahadrauf/Desktop/Research/medical_shape_display/data/strain_tests/"
file_name = "20240729_18_07_32_test_stepper_hv507_ad7124_onboardlogging_nodelaydx=1.5_speed=1s_V=300_f=10.npy"

time_loadcell, load_cell, stage, relay_output, idx_stop_moving, idx_relay_input_rise, idx_relay_input_fall = import_data(file_loc + file_name,
                                                                                                                         zero_load_cell=False,
                                                                                                                         convert_to_kPa=True,
                                                                                                                         width_pattern=2,
                                                                                                                         length_pattern=55.5)

fig = plt.figure(layout='constrained', figsize=(12, 4.5))
fig_top, fig_bot = fig.subfigures(nrows=1, ncols=2, width_ratios=[2, 1])
fig_bot2 = fig_bot

ax1_right = fig_top.add_subplot(111)
ax1 = fig_top.add_subplot(111, sharex=ax1_right, frameon=False)

color_left = "#1965B0"
color_right = "#DB060B"

ax1_right.plot(time_loadcell, relay_output, c=color_right, alpha=0.6)
ax1.plot(time_loadcell, load_cell, c=color_left)
ax1.set_xlabel("Time (s)")
ax1.set_xticks([0, 2, 4, 6])
ax1.set_ylabel("Shear Pressure (kPa)")
ax1.yaxis.label.set_color(color_left)
ax1.tick_params(axis='y', colors=color_left)
ax1_right.yaxis.tick_right()
ax1_right.yaxis.set_label_position("right")
ax1_right.set_ylabel("Voltage (V)")
ax1_right.yaxis.label.set_color(color_right)
ax1_right.tick_params(axis='y', colors=color_right)
ax1_right.grid(True)
ax1.grid(True)
ax1.set_yticks([-10, 0, 10, 20, 30])
ax1_right.grid(True, c=color_right, alpha=0.6)
ax1_right.set_yticks([0, 100, 200, 300])

ax2_right = fig_bot.add_subplot(211)
ax2 = fig_bot.add_subplot(211, sharex=ax2_right, frameon=False)
ax3_right = fig_bot2.add_subplot(212, sharey=ax2_right)
ax3 = fig_bot2.add_subplot(212, sharex=ax3_right, frameon=False)

# Fig. 5d
time_loadcell -= np.mean([time_loadcell[idx_relay_input_rise], time_loadcell[idx_relay_input_rise - 1]])
sos = butter(2, 250, btype='lowpass', fs=9300, output='sos')
load_cell_filtered = sosfiltfilt(sos, load_cell)
didx_plot = int(9300 * 0.0025)
line_active, = ax2.plot(1000*time_loadcell[idx_relay_input_rise - didx_plot:idx_relay_input_rise + didx_plot],
                        load_cell_filtered[idx_relay_input_rise - didx_plot:idx_relay_input_rise + didx_plot], c=color_left,
                        alpha=0.9)
ax2_right.plot(1000*time_loadcell[idx_relay_input_rise - didx_plot:idx_relay_input_rise + didx_plot],
               relay_output[idx_relay_input_rise - didx_plot:idx_relay_input_rise + didx_plot], c=color_right)

didx_fit_start = int(9300*0.0005)
engagement_times = []
didx_fit_ends = np.arange(int(9300 * 0.002), int(9300 * 0.02))
didx_fit_ends_valid = []
offset_when_voltage_enabled = np.nanmean(load_cell[idx_relay_input_rise - int(9300*0.003):idx_relay_input_rise])
for didx_fit_end in didx_fit_ends:
    time_fit = time_loadcell[idx_relay_input_rise + didx_fit_start:idx_relay_input_rise + didx_fit_end]
    slope, intercept, r, p, se = linregress(time_fit, load_cell[idx_relay_input_rise + didx_fit_start:idx_relay_input_rise + didx_fit_end])
    t_engage = (offset_when_voltage_enabled - intercept)/slope
    if t_engage > 0 and r*r > 0.8:
        engagement_times.append(t_engage)
        didx_fit_ends_valid.append(didx_fit_end)
t_engage = np.min(engagement_times)
didx_fit_end = didx_fit_ends_valid[np.argmin(engagement_times)]
time_fit = time_loadcell[idx_relay_input_rise + didx_fit_start:idx_relay_input_rise + didx_fit_end]
slope, intercept, r, p, se = linregress(time_fit, load_cell[idx_relay_input_rise + didx_fit_start:idx_relay_input_rise + didx_fit_end])
print("Engagement Time:", t_engage, t_engage - time_loadcell[idx_relay_input_rise], "R2", r*r)
print("Fitted between", time_loadcell[idx_relay_input_rise + didx_fit_start], "and", time_loadcell[idx_relay_input_rise + didx_fit_end])
time_fit = np.linspace(t_engage, min(time_fit[-1], time_loadcell[idx_relay_input_rise + didx_plot]))
ax2.axhline(offset_when_voltage_enabled, ls='--', c='k')
ax2.plot(time_fit*1e3, slope*time_fit + intercept, 'k', ls='--')
ax2.annotate(r"$t_{engage}$", xy=(t_engage*1e3 + 0.1, offset_when_voltage_enabled + 0.001), xytext=(t_engage*1e3 + 0.8, offset_when_voltage_enabled + 0.02), ha='left', va='center',
             fontsize=14, arrowprops=dict(arrowstyle="-|>", facecolor='black', relpos=(0, 0.5)))

ax2.grid(True)
ax2.set_xlabel("Time (ms) (+3 sec)")
ax2.set_ylabel("Shear Pressure (kPa)", y=0.35)
fig_bot.align_ylabels()
ax2.yaxis.label.set_color(color_left)
ax2.tick_params(axis='y', colors=color_left)
ax2_right.yaxis.tick_right()
ax2_right.yaxis.set_label_position("right")
ax2_right.set_ylabel("Voltage (V)")
ax2_right.yaxis.label.set_color(color_right)
ax2_right.tick_params(axis='y', colors=color_right)
ax2_right.grid(True, c=color_right, alpha=0.6)
ax2_right.grid(True, c=color_right, alpha=0.6)
ax2_right.set_yticks([0, 100, 200, 300])

# Plot fall time
print("Fell around t =", time_loadcell[idx_relay_input_fall])
didx_plot = int(9300*0.015)
ax3.plot(1e3*(time_loadcell[idx_relay_input_fall - didx_plot:idx_relay_input_fall + didx_plot] - time_loadcell[idx_relay_input_fall]),
         load_cell[idx_relay_input_fall - didx_plot:idx_relay_input_fall + didx_plot])

# The load cell plot is based on the original file, but I'm plotting just the relay voltage output based on another datasheet that recorded the voltage of the second electrode instead (to better show the release time)
file_name = "20240729_19_34_40_test_stepper_hv507_ad7124_onboardlogging_nodelaydx=1.5_speed=1s_V=300_f=1.npy"
time_loadcell, load_cell, stage, relay_output, idx_stop_moving, idx_relay_input_rise, idx_relay_input_fall = import_data(file_loc + file_name, zero_load_cell=False, convert_to_kPa=True, width_pattern=2, length_pattern=55.5)
ax3_right.plot(1e3*(time_loadcell[idx_relay_input_fall - didx_plot:idx_relay_input_fall + didx_plot] - time_loadcell[idx_relay_input_fall]),
               relay_output[idx_relay_input_fall - didx_plot:idx_relay_input_fall + didx_plot], c=color_right)
threshold = np.mean(load_cell[idx_relay_input_fall + didx_plot - 100:idx_relay_input_fall + didx_plot]) + \
            0.1*(np.max(load_cell[idx_relay_input_fall - 10:]) - np.mean(load_cell[idx_relay_input_fall + didx_plot - 100:idx_relay_input_fall + didx_plot]))
idx_cross_threshold = np.where((load_cell < threshold) & (time_loadcell > time_loadcell[idx_relay_input_fall]))[0][0]
time_cross_threshold = 1e3*(time_loadcell[idx_cross_threshold] - time_loadcell[idx_relay_input_fall])
print("Release time:", 1e3*(time_loadcell[idx_cross_threshold] - time_loadcell[idx_relay_input_fall]))
ymin, ymax = ax3.get_ylim()
print("Percentage", (ymin + threshold)/(ymax - ymin))
ax3.axvline(0, ls='--', c='k', ymax=(threshold - ymin)/(ymax - ymin))
ax3.axvline(time_cross_threshold, ls='--', c='k', ymax=(threshold - ymin)/(ymax - ymin))
ax3.axhline(threshold, ls='--', c='k')
ax3.annotate(r"$t_{release}$", xy=(0, (ymin + threshold)/2), xytext=(-4, (ymin + threshold)/2), ha='right', va='center',
             fontsize=14, arrowprops=dict(arrowstyle="-|>", facecolor='black'))
ax3.annotate("", xy=(time_cross_threshold, (ymin + threshold)/2), xytext=(time_cross_threshold + 4, (ymin + threshold)/2), ha='left', va='center',
             fontsize=14, arrowprops=dict(arrowstyle="-|>", facecolor='black'))

ax3.grid(True)
ax3.set_xlabel("Time (ms) (+4.5 sec)")
ax3.set_ylabel("Shear Pressure (kPa)", y=0.35)
ax3.yaxis.label.set_color(color_left)
ax3.tick_params(axis='y', colors=color_left)
ax3_right.yaxis.tick_right()
ax3_right.yaxis.set_label_position("right")
ax3_right.set_ylabel("Voltage (V)")
ax3_right.yaxis.label.set_color(color_right)
ax3_right.tick_params(axis='y', colors=color_right)
ax3_right.grid(True, c=color_right, alpha=0.6)
ax3_right.set_yticks([0, 100, 200, 300])

t1 = plt.text(0, 1, "(c)", fontsize=16, transform=fig.transFigure, horizontalalignment='left', verticalalignment='top')
t2 = plt.text(0.655, 1, "(d)", fontsize=16, transform=fig.transFigure, horizontalalignment='left', verticalalignment='top')
t3 = plt.text(0.655, 0.5, "(e)", fontsize=16, transform=fig.transFigure, horizontalalignment='left', verticalalignment='top')
t1.set_in_layout(False)
t2.set_in_layout(False)
t3.set_in_layout(False)

# plt.savefig(save_folder + timestamp + ".png", dpi=300)
# # plt.savefig("figures/" + timestamp + ".svg")
# plt.savefig(save_folder + timestamp + ".pdf")
plt.show()
