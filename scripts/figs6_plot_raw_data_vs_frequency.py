import numpy as np
import matplotlib.pyplot as plt
from import_strain_test_data import import_data
from calculate_timing_stats import calculate_engagement_time, calculate_release_time, calculate_Fpreload
from scipy.signal import butter, sosfiltfilt
from datetime import datetime

now = datetime.now()
name_clarifier = "_plot_strain_tests_raw_V=300_vs_freq"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
save_folder = "C:/Users/ahadrauf/Desktop/Research/latex/electroadhesive_dynamics_paper/figures/"

plt.rcParams["font.family"] = "Arial"
plt.rc('font', size=16)  # controls default text size
plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=16)  # fontsize of the y tick labels
plt.rc('legend', fontsize=16)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

file_loc = "../data/"
file_names = ["20241029_17_47_42_bts_dx=1.5_v=0.5_Fpre=-1g_V=300_f=0.1_w=2_brass",
              "20241029_17_53_39_bts_dx=1.5_v=0.5_Fpre=-1g_V=300_f=10_w=2_brass",
              "20241029_17_25_32_bts_dx=1.5_v=0.5_Fpre=-1g_V=300_f=100_w=2_brass",
              "20241027_18_31_57_bts_dx=1.5_v=0.5_Fpre=-1g_V=300_f=1000_w=2_brass",
              "20241029_17_33_51_bts_dx=1.5_v=0.5_Fpre=-1g_V=300_f=5000_w=2_brass",
              ]

w_pin = 2
depth_pin = 2

nx, ny = 5, 3
fig, axs = plt.subplots(nx, ny, layout='constrained', figsize=(10, 11))
frequencies = ["DC", "10 Hz", "100 Hz", "1000 Hz", "10000 Hz"]
for i, file_name in enumerate(file_names):
    ax1_right, ax2_right, ax3_right = axs[i]
    time_loadcell, load_cell, relay_input, relay_output, idx_stop_moving, idx_relay_input_rise, idx_relay_input_fall = import_data(file_loc + file_name + ".npy",
                                                                                                                                   width_pattern=w_pin,
                                                                                                                                   zero_load_cell=False,
                                                                                                                                   convert_to_kPa=False)
    sos = butter(5, (48, 52), btype='bandstop', fs=12800, output='sos')
    load_cell_filtered = sosfiltfilt(sos, load_cell)
    sos = butter(3, 250, btype='lowpass', fs=12800, output='sos')
    load_cell_filtered = sosfiltfilt(sos, load_cell_filtered)

    F_preload = calculate_Fpreload(file_loc + file_name + ".npy", w_pin=2e-3, depth_pin=2e-3)
    print("F_preload", F_preload)
    engagement_time = calculate_engagement_time(file_loc + file_name + ".npy")
    release_time = calculate_release_time(file_loc + file_name + ".npy")
    if i == 0:
        ax1 = fig.add_subplot(nx, ny, 3*i + 1, sharex=ax1_right, frameon=False)
        ax2 = fig.add_subplot(nx, ny, 3*i + 2, sharex=ax2_right, frameon=False)
        ax3 = fig.add_subplot(nx, ny, 3*i + 3, sharex=ax3_right, frameon=False)
        top_ax1, top_ax2, top_ax3 = ax1, ax2, ax3
    else:
        ax1 = fig.add_subplot(nx, ny, 3*i + 1, sharex=ax1_right, sharey=top_ax1, frameon=False)
        ax2 = fig.add_subplot(nx, ny, 3*i + 2, sharex=ax2_right, sharey=top_ax2, frameon=False)
        ax3 = fig.add_subplot(nx, ny, 3*i + 3, sharex=ax3_right, sharey=top_ax3, frameon=False)

    ax1.plot(time_loadcell, load_cell, color='tab:blue', zorder=89)
    ax1.plot(time_loadcell, load_cell_filtered, color='tab:orange', zorder=90)
    ax1_right.plot(time_loadcell, relay_output, color='tab:red', alpha=0.8)
    ax1.text(0.03, 0.03, frequencies[i], transform=ax1.transAxes, horizontalalignment='left', verticalalignment='bottom',
             bbox=dict(facecolor='white', alpha=0.75, linewidth=0), fontsize=14)
    ax2.text(0.03, 0.03, frequencies[i], transform=ax2.transAxes, horizontalalignment='left', verticalalignment='bottom',
             bbox=dict(facecolor='white', alpha=0.75, linewidth=0), fontsize=14)
    ax3.text(0.03, 0.03, frequencies[i], transform=ax3.transAxes, horizontalalignment='left', verticalalignment='bottom',
             bbox=dict(facecolor='white', alpha=0.75, linewidth=0), fontsize=14)

    didx_plot = int(12800*0.020)
    time_plot = 1000*(time_loadcell[idx_relay_input_rise - didx_plot:idx_relay_input_rise + didx_plot] - time_loadcell[idx_relay_input_rise])
    ax2.plot(time_plot,
             load_cell[idx_relay_input_rise - didx_plot:idx_relay_input_rise + didx_plot], color='tab:blue', zorder=89)
    ax2.plot(time_plot,
             load_cell_filtered[idx_relay_input_rise - didx_plot:idx_relay_input_rise + didx_plot], color='tab:orange', zorder=90)
    ax2_right.plot(time_plot,
                   relay_output[idx_relay_input_rise - didx_plot:idx_relay_input_rise + didx_plot], color='tab:red', alpha=0.8)

    didx_plot = int(12800*0.025)
    time_plot = 1000*(time_loadcell[idx_relay_input_fall - didx_plot:idx_relay_input_fall + didx_plot] - time_loadcell[idx_relay_input_fall])
    ax3.plot(time_plot,
             load_cell[idx_relay_input_fall - didx_plot:idx_relay_input_fall + didx_plot], color='tab:blue', zorder=89)
    ax3_right.plot(time_plot,
                   relay_output[idx_relay_input_fall - didx_plot:idx_relay_input_fall + didx_plot], color='tab:red', alpha=0.8)

    offset_when_voltage_enabled = np.nanmean(load_cell_filtered[idx_relay_input_rise - int(12800*0.01):idx_relay_input_rise - int(12800*0.00)])
    ax2.axhline(offset_when_voltage_enabled, ls='--', c='k', zorder=94)

    ax1_right.yaxis.tick_right()
    ax2_right.yaxis.tick_right()
    ax2_right.set_yticks([0, 100, 200, 300])
    ax3_right.yaxis.tick_right()
    ax1_right.yaxis.set_label_position("right")
    ax2_right.yaxis.set_label_position("right")
    ax3_right.yaxis.set_label_position("right")
    ax1_right.yaxis.label.set_color('tab:red')
    ax2_right.yaxis.label.set_color('tab:red')
    ax3_right.yaxis.label.set_color('tab:red')
    ax1_right.tick_params(axis='y', colors='tab:red')
    ax2_right.tick_params(axis='y', colors='tab:red')
    ax3_right.tick_params(axis='y', colors='tab:red')
    ax1.grid(True)
    ax2.grid(True)
    ax3.grid(True)

    if i == 0:
        ax1.set_title("Entire Run")
        ax2.set_title("Engagement Time")
        ax3.set_title("Release Time")
    if i == len(file_names) - 1:
        ax1.set_xlabel("Time (s)")
        ax2.set_xlabel("Time (ms + 3 sec)")
        ax3.set_xlabel("Time (ms + 4.5 sec)")
    if i == ny//2 + 1:
        ax1.set_ylabel("Load Force (N)", fontsize=18)
        ax3_right.set_ylabel("Voltage (V)", fontsize=18)

# plt.savefig(save_folder + timestamp + ".png", dpi=300)
# # # plt.savefig("figures/" + timestamp + ".svg")
# plt.savefig(save_folder + timestamp + ".pdf")

plt.show()
