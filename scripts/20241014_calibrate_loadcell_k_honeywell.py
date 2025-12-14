import numpy as np
import matplotlib.pyplot as plt
from import_strain_test_data_20240726 import import_data
from scipy.signal import savgol_filter, butter, sosfilt
from scipy.optimize import minimize
from sklearn.metrics import r2_score
from general.general_dsp import fft, power
from scipy.stats import linregress

file_loc = "../data/strain_tests/"

# file_names = {0.1: "20240814_21_24_27_calibrate_loadcell_k_dx=1_speed=0.1mms",
#               0.2: "20240814_21_29_47_calibrate_loadcell_k_dx=1_speed=0.2mms",
#               0.3: "20240814_21_32_07_calibrate_loadcell_k_dx=1_speed=0.3mms",
#               0.4: "20240814_21_34_37_calibrate_loadcell_k_dx=1_speed=0.4mms",
#               0.5: "20240814_21_20_19_calibrate_loadcell_k_dx=1_speed=0.5mms",
#               0.6: "20240814_21_21_14_calibrate_loadcell_k_dx=1_speed=0.6mms",
#               0.7: "20240814_21_36_26_calibrate_loadcell_k_dx=1_speed=0.7mms",
#               0.8: "20240814_21_22_28_calibrate_loadcell_k_dx=1_speed=0.8mms",
#               0.9: "20240814_21_23_03_calibrate_loadcell_k_dx=1_speed=0.9mms",
#               1: "20240814_21_18_17_calibrate_loadcell_k_dx=1_speed=1mms",
#               2: "20240814_21_38_56_calibrate_loadcell_k_dx=1_speed=2mms",
#               3: "20240814_21_39_18_calibrate_loadcell_k_dx=1_speed=3mms",
#               4: "20240814_21_39_36_calibrate_loadcell_k_dx=1_speed=4mms",
#               5: "20240814_21_39_56_calibrate_loadcell_k_dx=1_speed=5mms",
#               6: "20240814_21_24_06_calibrate_loadcell_k_dx=1_speed=6mms",
#               7: "20240814_21_40_14_calibrate_loadcell_k_dx=1_speed=7mms",
#               8: "20240814_21_40_32_calibrate_loadcell_k_dx=1_speed=8mms"}
file_names = {0.1: ["20241010_18_14_29_calibrate_loadcell_k_dx=1_speed=0.1mms",
                    "20241010_18_12_44_calibrate_loadcell_k_dx=1_speed=0.1mms",
                    "20241010_18_10_58_calibrate_loadcell_k_dx=1_speed=0.1mms",
                    "20241010_18_09_12_calibrate_loadcell_k_dx=1_speed=0.1mms"],
              0.5: ["20241010_18_07_22_calibrate_loadcell_k_dx=1_speed=0.5mms",
                    "20241010_18_07_01_calibrate_loadcell_k_dx=1_speed=0.5mms",
                    "20241010_18_06_40_calibrate_loadcell_k_dx=1_speed=0.5mms",
                    "20241010_18_06_19_calibrate_loadcell_k_dx=1_speed=0.5mms"],
              1: ["20241010_18_08_42_calibrate_loadcell_k_dx=1_speed=1mms",
                  "20241010_18_08_21_calibrate_loadcell_k_dx=1_speed=1mms",
                  "20241010_18_08_31_calibrate_loadcell_k_dx=1_speed=1mms",
                  "20241010_18_08_52_calibrate_loadcell_k_dx=1_speed=1mms"]}

all_speeds = sorted(file_names.keys())
all_k_mean = []
all_k_std = []
all_zero_std = []
all_max_noise = []
for speed in sorted(file_names.keys()):
    noise_freqs = []
    spring_constants = []
    zero_stds = []
    for file_name in file_names[speed]:
        data = np.load(file_loc + file_name + ".npy")
        print(np.shape(data))
        time_loadcell, load_cell, stage = data[:, 0], data[:, 1], data[:, 2]

        # Process data to remove noise/anomalies
        time_loadcell -= time_loadcell[0]

        # Delete NaN values
        to_delete = []
        ks, ms, bs = [], [], []
        for idx in range(len(time_loadcell)):
            if np.isnan(time_loadcell[idx]) or np.isnan(load_cell[idx]) or load_cell[idx] > 50 or time_loadcell[idx] > 1e8:
                to_delete.append(idx)
            if idx >= 1 and abs(time_loadcell[idx] - time_loadcell[idx - 1]) >= 1000:
                to_delete.append(idx)
            elif idx >= 1 and time_loadcell[idx] <= time_loadcell[idx - 1]:
                to_delete.append(idx)
            elif load_cell[idx] <= -10 or load_cell[idx] > 20:
                to_delete.append(idx)
        print("Deleted {} points".format(len(to_delete)))
        time_loadcell = np.delete(time_loadcell, to_delete)
        load_cell = np.delete(load_cell, to_delete)
        stage = np.delete(stage, to_delete)
        time_loadcell /= 1e6

        zero_mean = np.nanmean(load_cell[:100])
        zero_std = np.nanstd(load_cell[:100])
        idx_max_fit = np.nanargmax(load_cell)
        idx_fit = np.where(load_cell[:idx_max_fit] <= zero_mean + 3 * zero_std)[0][-1] + 1
        zero_mean = np.nanmean(load_cell[:idx_fit // 2])
        zero_std = np.nanstd(load_cell[:idx_fit // 2])
        zero_fft_freq, zero_fft = fft(time_loadcell[:idx_fit // 2], load_cell[:idx_fit // 2], oversample=5)

        sol = linregress(time_loadcell[idx_fit:idx_max_fit], load_cell[idx_fit:idx_max_fit])
        slope = abs(sol.slope)
        r_value = sol.rvalue
        spring_constants.append(slope / speed)

        zero_stds.append(zero_std)
        zero_fft_freq, zero_fft = zero_fft_freq[:len(zero_fft) // 2], zero_fft[:len(zero_fft) // 2]
        noise_freqs.append(zero_fft_freq[np.argmax(power(zero_fft)[25:]) + 25])
        disp_string = 'Zero Mean: {:0.3f}, Zero Std: {:0.3f}, Fit at t = {:0.3f}, Max Noise at: {:0.3f} Hz --> k = {:0.3f}, r2 = {:0.3f}' + \
            " (({:0.3f}, {:0.3f}) --> ({:0.3f}, {:0.3f}))"
        print(disp_string.format(zero_mean, zero_std, time_loadcell[idx_fit], noise_freqs[-1], spring_constants[-1], r_value**2,
                                 time_loadcell[idx_fit], load_cell[idx_fit], time_loadcell[idx_max_fit], load_cell[idx_max_fit]))
        # print("Zero Mean:", zero_mean, "Zero Std:", zero_std, "Time Fit:", time_fit[idx_fit], "Max Noise At:", noise_freqs[-1], 'Hz')

        # print(np.max(power(zero_fft)))
        # plt.plot(zero_fft_freq, power(zero_fft))
        plt.plot(time_loadcell * speed, load_cell)

        print("k for speed =", speed, ":", np.nanmean(spring_constants), np.nanstd(spring_constants), spring_constants)
        all_k_mean.append(np.nanmean(spring_constants))
        all_k_std.append(np.nanstd(spring_constants))
        all_zero_std.append(np.nanmean(zero_stds))
    print("Average spring constant for speed =", speed, "=", np.mean(spring_constants), np.std(spring_constants))
    # plt.plot(time_loadcell, load_cell)
    # plt.plot(time_loadcell, stage)

print("All spring constants", np.mean(all_k_mean), all_k_mean, all_k_std)
# plt.errorbar(all_speeds, all_k_mean, all_k_std, capsize=3)
plt.grid(True)
# plt.legend(["Real", "Fitted"])
# plt.xlabel("Time (s)")
# plt.ylabel("Load Cell (N)")
plt.show()
