import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import matplotlib.patches as mpatches

now = datetime.now()
name_clarifier = "_adc_timing_oscilloscope"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
save_folder = "../figures/"
start_time = datetime.now()

plt.rcParams["font.family"] = "Arial"
plt.rc('font', size=16)  # controls default text size
plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=16)  # fontsize of the y tick labels
plt.rc('legend', fontsize=16)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
colors = ["#1964B0", "#F1932D", "#4DB264", "#DB060B"]

folder_loc = "../data/"
file_name_dataready = "20250416_oscilloscope_dataready.csv"
file_name_motorstep = "20250416_oscilloscope_motorstep.csv"

fig, axs = plt.subplots(1, 2, layout='constrained', figsize=(12, 5))
ax1_right, ax2_right = axs
ax1 = fig.add_subplot(121, sharex=ax1_right, frameon=False)
ax2 = fig.add_subplot(122, sharex=ax2_right, frameon=False)

data = np.genfromtxt(folder_loc + file_name_dataready, delimiter=',')
time, vin, dataready, step = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
start_idx = np.where(vin > 3.3 * 0.75)[0][0]
time -= time[start_idx]

# Sync motor timing signals with microcontroller timing signals
data2 = np.genfromtxt(folder_loc + file_name_motorstep, delimiter=',')
time2, motor2, step2 = data2[:, 0], data2[:, 1], data2[:, 3]
cutoff_idx = 50
time2, motor2, step2 = time2[cutoff_idx:], motor2[cutoff_idx:], step2[cutoff_idx:]
step2_rise_idx = np.where(step2 > 1.5)[0][0]  # from datasheet, Table 6.5: https://www.ti.com/lit/ds/symlink/drv8711.pdf
time2_shift = time2[step2_rise_idx]

step1_skip_idx = np.where(step > 3.3 * 0.8)[0][0]
step1_rise_idx = np.where((time > time[step1_skip_idx]) & (step < 1.5))[0][0]  # from datasheet, Table 6.5: https://www.ti.com/lit/ds/symlink/drv8711.pdf
time1_shift = time[step1_rise_idx]
time2 = time2 - time2_shift + time1_shift
print("Time 1 shift", time1_shift)
print("Time 2 shift", time2_shift)

min_time2_idx = np.where(time2 > time[0])[0][0]
max_time2_idx = np.where(time2 > time[-1])[0][0] - 1

plot_hvinterrupt, = ax1.plot(1e6 * time, vin, c=colors[0], label="High Voltage Interrupt")
plot_dataready, = ax1.plot(1e6 * time, dataready, c=colors[2], label="ADC Data Ready")
plot_motorinput, = ax1_right.plot(1e6 * time2[min_time2_idx:max_time2_idx], motor2[min_time2_idx:max_time2_idx], c=colors[3], label="Motor Input Voltage")
ax1.set_xlabel("Time (μs)")
ax1.set_ylabel("Logic Voltages (V)")
ax1.grid(True)
ax1.legend([plot_hvinterrupt, plot_dataready, plot_motorinput], ["High Voltage Interrupt", "ADC Data Ready", "Motor Input Voltage"], loc='lower right')

ax1_right.yaxis.tick_right()
ax1_right.yaxis.set_label_position("right")
ax1_right.set_ylabel("Motor Input Voltage (V)")
ax1_right.yaxis.label.set_color('tab:red')
ax1_right.tick_params(axis='y', colors='tab:red')

ax1.set_ylim([-0.5, 5.5])
ax1_right.set_ylim([-14.5 * 0.5 / 5.5, 14.5])

# Plot zoomed in figure
file_name_dataready = "20250416_oscilloscope_dataready_zoomedin.csv"
file_name_motorstep = "20250416_oscilloscope_motorstep_zoomedin.csv"
data = np.genfromtxt(folder_loc + file_name_dataready, delimiter=',')
time, vin, dataready, step = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
start_idx = np.where(vin > 3.3 * 0.95)[0][0]
time -= time[start_idx]

# Sync motor timing signals with microcontroller timing signals
data2 = np.genfromtxt(folder_loc + file_name_motorstep, delimiter=',')
time2, motor2, step2 = data2[:, 0], data2[:, 1], data2[:, 3]
step2_rise_idx = np.where(step2 > 1.5)[0][0]  # from datasheet, Table 6.5: https://www.ti.com/lit/ds/symlink/drv8711.pdf
time2_shift = time2[step2_rise_idx]

step1_rise_idx = np.where(step < 3.3 * 0.5)[0][0]
time1_shift = time[step1_rise_idx]
time2 = time2 - time2_shift + time1_shift
print("Time 1 shift", time1_shift)
print("Time 2 shift", time2_shift)

max_time2_idx = np.where(time2 > time[-1])[0][0] - 1
avg_base_motor2 = np.mean(motor2[-100:])
time2 = np.append([min(time)], time2)
motor2 = np.append([avg_base_motor2], motor2)

start_idx = np.where(vin > 3.3 * 0.9)[0][0]
time -= time[start_idx]
end_idx = np.where(dataready > 3.3 * 0.9)[0][0]
delay = 1e6 * (time[end_idx] - time[start_idx])
print("Delay = ", time[end_idx] - time[start_idx])
arrow = mpatches.FancyArrowPatch((time[start_idx] * 1e6, 3.3 * 0.8), (time[end_idx] * 1e6, 3.3 * 0.8), arrowstyle='<->,head_width=.15', mutation_scale=20)
ax2.add_patch(arrow)
ax2.annotate("Delay = {:0.2f} μs".format(delay), (.5, 0), xycoords=arrow, ha='center', va='top')
plot_hvinterrupt, = ax2.plot(1e6 * time, vin, c=colors[0], label="High Voltage Interrupt")
plot_dataready, = ax2.plot(1e6 * time, dataready, c=colors[2], label="ADC Data Ready")
plot_motorinput, = ax2_right.plot(1e6 * time2[:max_time2_idx], motor2[:max_time2_idx], c=colors[3], label="Motor Input Voltage")
ax2.set_xlabel("Time (μs)")
ax2.set_ylabel("Logic Voltages (V)")
ax2.grid(True)
ax2.legend([plot_hvinterrupt, plot_dataready, plot_motorinput], ["High Voltage Interrupt", "ADC Data Ready", "Motor Input Voltage"], loc='lower right')

ax2_right.yaxis.tick_right()
ax2_right.yaxis.set_label_position("right")
ax2_right.set_ylabel("Motor Input Voltage (V)")
ax2_right.yaxis.label.set_color('tab:red')
ax2_right.tick_params(axis='y', colors='tab:red')

ax2.set_ylim([-0.5, 5.5])
ax2_right.set_ylim([-14.5 * 0.5 / 5.5, 14.5])

fig.text(0.00, .99, "(a)", transform=fig.transFigure, horizontalalignment='left', verticalalignment='top')
fig.text(0.5, .99, "(b)", transform=fig.transFigure, horizontalalignment='left', verticalalignment='top')

# fig.savefig(save_folder + timestamp + ".png", dpi=300)
# # # plt.savefig("figures/" + timestamp + ".svg")
# fig.savefig(save_folder + timestamp + ".pdf")

plt.show()
