import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from import_strain_test_data import import_data
from scipy.signal import butter, sosfiltfilt

now = datetime.now()
name_clarifier = "_noise_fft"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
save_folder = "../data/"
start_time = datetime.now()

plt.rcParams["font.family"] = "Arial"
plt.rc('font', size=16)  # controls default text size
plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=16)  # fontsize of the y tick labels
plt.rc('legend', fontsize=16)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

fig, axs = plt.subplots(1, 2, layout='constrained', figsize=(12, 5))
ax1, ax2 = axs

file_loc = "../data/"
file_names = ["20241107_14_50_58_bts_dx=1.5_v=0.5_Fpre=-1g_V=300_f=0.1_w=2_brass",
              "20241107_14_49_27_bts_dx=1.5_v=0.5_Fpre=-1g_V=300_f=0.1_w=2_brass",
              "20241107_14_48_44_bts_dx=1.5_v=0.5_Fpre=-1g_V=300_f=0.1_w=2_brass",
              "20241107_14_48_01_bts_dx=1.5_v=0.5_Fpre=-1g_V=300_f=0.1_w=2_brass",
              "20241107_14_47_17_bts_dx=1.5_v=0.5_Fpre=-1g_V=300_f=0.1_w=2_brass",
              "20241107_14_45_51_bts_dx=1.5_v=0.5_Fpre=-1g_V=300_f=0.1_w=2_brass",
              "20241107_15_06_58_bts_dx=1.5_v=0.5_Fpre=-1g_V=200_f=10_w=2_brass",
              "20241107_15_06_15_bts_dx=1.5_v=0.5_Fpre=-1g_V=200_f=10_w=2_brass",
              "20241107_15_05_05_bts_dx=1.5_v=0.5_Fpre=-1g_V=200_f=10_w=2_brass",
              "20241107_15_04_17_bts_dx=1.5_v=0.5_Fpre=-1g_V=200_f=10_w=2_brass",
              "20241107_15_03_28_bts_dx=1.5_v=0.5_Fpre=-1g_V=200_f=10_w=2_brass",
              "20241107_15_02_45_bts_dx=1.5_v=0.5_Fpre=-1g_V=200_f=10_w=2_brass",
              "20241107_15_00_35_bts_dx=1.5_v=0.5_Fpre=-1g_V=200_f=0.1_w=2_brass",
              "20241107_14_59_52_bts_dx=1.5_v=0.5_Fpre=-1g_V=200_f=0.1_w=2_brass",
              "20241107_14_59_08_bts_dx=1.5_v=0.5_Fpre=-1g_V=200_f=0.1_w=2_brass",
              "20241107_14_58_07_bts_dx=1.5_v=0.5_Fpre=-1g_V=200_f=0.1_w=2_brass",
              "20241107_14_57_16_bts_dx=1.5_v=0.5_Fpre=-1g_V=200_f=0.1_w=2_brass",
              "20241107_14_56_33_bts_dx=1.5_v=0.5_Fpre=-1g_V=200_f=0.1_w=2_brass",
              "20241107_15_08_21_bts_dx=1.5_v=0.5_Fpre=-1g_V=200_f=10000_w=2_brass",
              "20241107_15_09_04_bts_dx=1.5_v=0.5_Fpre=-1g_V=200_f=10000_w=2_brass"]
print("# Files:", len(file_names))


def fft(t, x, oversample=1):
    avg_dt = np.nanmean(t[1:] - t[:-1])
    return np.abs(np.fft.fftfreq(np.size(t) * oversample, d=avg_dt)), np.fft.fft(x, n=np.size(t) * oversample)


def power(x):
    return np.square(np.abs(x)) / np.size(x)

all_k_mean = []
all_k_std = []
all_zero_std = []
all_max_noise = []
fft_powers = []
fft_powers_after_filter = []
filtered_times = np.linspace(1, 2, 12800 * 2)
filtered_load_cell_readings = []
for idx, file_name in enumerate(file_names):
    w_pin = 2
    depth_pin = 2
    time_loadcell, load_cell, relay_input, relay_output, idx_stop_moving, idx_relay_input_rise, idx_relay_input_fall = import_data(file_loc + file_name + ".npy",
                                                                                                                                   zero_load_cell=False,
                                                                                                                                   convert_to_kPa=False)

    noise_freqs = []
    spring_constants = []
    zero_stds = []
    oversample = 1

    idx_relay_input_rise = max(np.where(relay_output > 50)[0][0], np.where(relay_input > 0.5)[0][0])
    offset_when_voltage_enabled = np.nanmean(load_cell[idx_relay_input_rise - int(12800 * 0.102):idx_relay_input_rise - int(12800 * 0.002)])
    min_fit_idx = int(12800 * 1)
    time_loadcell_fit = time_loadcell[min_fit_idx:idx_relay_input_rise - int(12800 * 0.002)]

    loadcell_fit = load_cell[min_fit_idx:idx_relay_input_rise - int(12800 * 0.002)]
    zero_fft_freq, zero_fft = fft(time_loadcell_fit, loadcell_fit, oversample=oversample)
    zero_fft_freq, zero_fft = zero_fft_freq[:len(zero_fft) // 2], zero_fft[:len(zero_fft) // 2]
    noise_freqs.append(zero_fft_freq[np.argmax(power(zero_fft)[25:]) + 25])
    all_zero_std.append(np.std(loadcell_fit))
    line_timedata, = ax1.plot(time_loadcell_fit, loadcell_fit)

    idx_max_freq = np.where(zero_fft_freq > 1000)[0][0]
    fft_powers.append(power(zero_fft[1:idx_max_freq]))

    sos = butter(3, 250, btype='lowpass', fs=12800, output='sos')
    load_cell_filtered = sosfiltfilt(sos, load_cell)
    loadcell_fit = load_cell_filtered[min_fit_idx:idx_relay_input_rise - int(12800 * 0.002)]
    zero_fft_freq, zero_fft = fft(time_loadcell_fit, loadcell_fit, oversample=oversample)
    zero_fft_freq, zero_fft = zero_fft_freq[:len(zero_fft) // 2], zero_fft[:len(zero_fft) // 2]
    fft_powers_after_filter.append(power(zero_fft[1:idx_max_freq]))

    filtered_load_cell_readings.append(np.interp(filtered_times, time_loadcell_fit, loadcell_fit))

    if idx == len(file_names) - 1:
        line_filtereddata, = ax1.plot(time_loadcell_fit, loadcell_fit, color='k')
        ax1.legend([line_timedata, line_filtereddata], ["Original Data", "Data After 250 Hz Low-Pass"])

plot_fft_before_lowpass, = ax2.plot(zero_fft_freq[1:idx_max_freq], np.mean(fft_powers, axis=0), color='tab:blue')
factor = np.divide(np.mean(fft_powers, axis=0) + np.std(fft_powers, axis=0), np.mean(fft_powers, axis=0))
ax2.fill_between(zero_fft_freq[1:idx_max_freq], np.mean(fft_powers, axis=0) / factor,
                 np.mean(fft_powers, axis=0) * factor, color='tab:blue', alpha=0.25)

factor = np.divide(np.mean(fft_powers_after_filter, axis=0) + np.std(fft_powers_after_filter, axis=0), np.mean(fft_powers_after_filter, axis=0))
plot_fft_after_lowpass, = ax2.plot(zero_fft_freq[1:idx_max_freq], np.mean(fft_powers_after_filter, axis=0), color='tab:orange')
ax2.fill_between(zero_fft_freq[1:idx_max_freq], np.mean(fft_powers_after_filter, axis=0) / factor,
                 np.mean(fft_powers_after_filter, axis=0) * factor, color='tab:orange', alpha=0.25)

print("Average standard deviation:", np.mean(all_zero_std), np.std(all_zero_std))
ax1.grid(True)
dtime = time_loadcell_fit[0]
ax1.set_xlabel("Time (s)")
ax1.set_ylabel("Load Cell Reading (N)")

ax2.axvline(250, ls='--', c='k')
ax2.set_xticks([0, 250, 500, 750, 1000])
ax2.grid(True)
ax2.set_xlabel("Frequency (Hz)")
ax2.set_ylabel("FFT Power")
ax2.set_yscale('log')
plt.legend([plot_fft_before_lowpass, plot_fft_after_lowpass], ["FFT of Original Data", "FFT After 250 Hz Low-Pass"])

fig.text(0.0, 0.99, "(a)", transform=fig.transFigure, horizontalalignment='left', verticalalignment='top')
fig.text(0.505, 0.99, "(b)", transform=fig.transFigure, horizontalalignment='left', verticalalignment='top')

# fig.savefig(save_folder + timestamp + ".png", dpi=300)
# # # plt.savefig("figures/" + timestamp + ".svg")
# fig.savefig(save_folder + timestamp + ".pdf")

plt.show()

