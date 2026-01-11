import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.signal import savgol_filter
from scipy.integrate import cumulative_trapezoid

now = datetime.now()
name_clarifier = "_current_vs_time_sample"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
plt.rcParams["font.family"] = "Arial"
plt.rc('font', size=20)  # controls default text size
plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=16)  # fontsize of the y tick labels
plt.rc('legend', fontsize=16)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

folder_loc = "C:/Users/ahadrauf/Desktop/Research/medical_shape_display/data/oscilloscope/20240826_oscilloscope_hv507_differential/"
file_name = "20240826_18_26_47_oscilloscope_hv507_diff_w=6.35_100.01kOhm_300V_1000Hz"
save_folder = "C:/Users/ahadrauf/Desktop/Research/latex/electroadhesive_dynamics_paper/figures/"


def process_data(file_loc):
    data = np.genfromtxt(file_loc + ".csv", delimiter=",")
    time, v1, v2, c1 = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    c1 /= 100.01e3  # --> mA

    print(np.shape(data))
    window_length = 21
    v1 = savgol_filter(v1, window_length, polyorder=2)
    v2 = savgol_filter(v2, window_length, polyorder=2)
    c1 = savgol_filter(c1, window_length, polyorder=2)

    idx_v1_rise1 = np.argmin(v1[:500])  # np.where(v1 > threshold_low)[0][0] - 1
    print(time[idx_v1_rise1])
    time -= time[idx_v1_rise1]
    idx_v1_fall1 = np.argmin(np.abs(time - 500e-6))  # for 1000 Hz signal, #np.where((v1 >= threshold_low) & (time > time[idx_v1_rise1] + 100))[0][0] - 1
    idx_v1_rise2 = np.argmin(np.abs(time - 1000e-6))  # np.where(v1 > threshold_high)[0][-1] + 1

    q1 = cumulative_trapezoid(c1, time, initial=0)
    c1_zero = (q1[idx_v1_rise2] - q1[idx_v1_rise1])/(time[idx_v1_rise2] - time[idx_v1_rise1])
    c1 -= c1_zero

    c1_first_avg = np.mean(c1[idx_v1_rise1 + 200:idx_v1_fall1 - 100])
    c1_second_avg = np.mean(c1[idx_v1_fall1 + 200:idx_v1_rise2 - 100])
    q1 = cumulative_trapezoid(c1, time, initial=0)
    q1 -= q1[idx_v1_rise1]
    q1_total = np.trapz(c1[idx_v1_rise1:idx_v1_rise2], time[idx_v1_rise1:idx_v1_rise2])
    print("Total charge integration:", q1_total, q1[idx_v1_rise2], time[idx_v1_rise2])
    dc1 = np.hstack([np.nan, np.divide(c1[2:] - c1[:-2], time[2:] - time[:-2]), np.nan])

    print("Time rise 1:", time[idx_v1_rise1], "fall 1", time[idx_v1_fall1], "Zero current:", c1_first_avg, "(b/w t =", time[idx_v1_rise1 + 200], time[idx_v1_fall1 - 100])
    print("Time rise 2:", time[idx_v1_rise2], "Zero current:", c1_second_avg, "(b/w t =", time[idx_v1_fall1 + 200], time[idx_v1_rise2 - 1000])

    power = np.multiply(v1[idx_v1_rise1:idx_v1_rise2], c1[idx_v1_rise1:idx_v1_rise2])
    avg_power = np.mean(power)
    print("Average power over cycle * 2", 2*avg_power*1000, "mW")

    return time, v1, v2, c1, q1


time, v1, v2, c1, q1 = process_data(folder_loc + file_name)
fig, ax1 = plt.subplots(layout='constrained', figsize=(6, 4))

colors = ["#1964B0", "#F1932D", "#4DB264", "#DB060B"]
line_v1, = ax1.plot(time*1e6, v1, c=colors[0], zorder=99)
line_v2, = ax1.plot(time*1e6, v2, c=colors[2])
ax1.set_xlabel("Time (µs)")
ax1.set_ylabel("Voltage (V)")
print("Max voltage", np.max(v1))

ax1_right = ax1.twinx()
line_c1, = ax1_right.plot(time*1e6, c1*1e3, colors[3])
ax1_right.set_ylabel("Current (mA)", color=colors[3])
ax1_right.tick_params(axis='y', labelcolor=colors[3])

ax1.grid(True)
ax1_right.grid(True, axis='y', color=colors[3], alpha=0.4)
plt.legend([line_v1, line_v2, line_c1], [r"$V_1$", r"$V_2$", r"$I_1$"], loc='lower right', fontsize=16)

# plt.savefig(save_folder + timestamp + ".png", dpi=300)
# # plt.savefig(save_folder + timestamp + ".svg")
# plt.savefig(save_folder + timestamp + ".pdf")
plt.show()
