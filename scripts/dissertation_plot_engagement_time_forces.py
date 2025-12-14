import numpy as np
import matplotlib.pyplot as plt
from sim_engagement_time import *
from datetime import datetime

now = datetime.now()
name_clarifier = "_plot_engagement_time_forces"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
save_folder = "../figures/"

plt.rcParams["font.family"] = "Arial"
plt.rc('font', size=16)  # controls default text size
plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=16)  # fontsize of the y tick labels
plt.rc('legend', fontsize=16)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
L_dielectric = 55.5e-3
w_pin = 2e-3
depth_pin = 2e-3
V_max = 300
k_dielectric = 54.2
t_dielectric = 24e-6
frequency = 1000
period = 1 / 2 / frequency
Fext = 0.125  # 0.04273761595304056
if w_pin < 1e-3:
    max_time_step = 0.1e-6
elif w_pin > 30e-3:
    max_time_step = 10e-6
else:
    max_time_step = 1e-6

T, Y, dY, X, dX, events, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic = sim_engagement_time(L_dielectric=L_dielectric, w_pin=w_pin,
                                                                                                       depth_pin=depth_pin, V_max=V_max,
                                                                                                       t_dielectric=t_dielectric, period=period,
                                                                                                       Fext=Fext, k_dielectric=k_dielectric,
                                                                                                       max_time_step=max_time_step)
print("Events:", events)
fig, axs = plt.subplots(1, 2, layout='constrained', figsize=(12, 4))
ax1, ax2 = axs
ax1_right = ax1.twinx()

# fig = plt.figure(layout='compressed', figsize=(12, 4))
# gs = plt.GridSpec(1, 2)
# ax1 = fig.add_subplot(gs[0])
# ax1_right = ax1.twinx()
# ax2 = fig.add_subplot(gs[1])

line_pos, = ax1.plot(T * 1e6, Y * 1e6, c='k', label=r"$T_{air}$")
# ax1.axhline(X[0] * 1e6, c='k', ls='--')
ax1.set_xlabel("Time (μs)")
ax1.set_ylabel("Position (μm)")
# ax.axhline(t_air_final * 1e6, ls='--', c='k')
line_vel, = ax1_right.plot(T * 1e6, dY * 1e3, c='r', label=r"$\dot{T}_{air}$")
# ax1.legend(loc='upper right')
# ax1_right.legend(loc='lower right')
ax1_right.legend([line_pos, line_vel], [r"$T_{air}$", r"$\dot{T}_{air}$"])
ax1_right.set_ylabel("Velocity (mm/s)")
ax1.grid(True)
ax1_right.yaxis.label.set_color('red')
ax1_right.tick_params(axis='y', colors='red')

colors = ['#0077BB', '#33BBEE', '#009988', '#EE7733', '#CC3311', '#EE3377']
colors = ['#4477AA', '#66CCEE', '#228833', '#CCBB44', '#EE6677', '#AA3377']
colors = ['#882E72', '#1965B0', '#7BAFDE', '#4EB265', '#EE8026', '#DC050C']
ax2.plot(T * 1e6, Fes, label=r'$F_{ea}$', color=colors[0])
ax2.plot(T * 1e6, Fk, label=r'$F_N$', color=colors[1])
ax2.plot(T * 1e6, Fb, label=r'$F_b$', color=colors[2])
ax2.plot(T * 1e6, Fconstant, label=r'$F_{preload}$', ls='--', color=colors[3])
# ax2.plot(T * 1e6, k_lc * X, label='Load Cell Spring')
# ax2.plot(T * 1e6, Flinear_static, label='Ffriction_s')
# ax2.plot(T * 1e6, Flinear_kinetic, label='Ffriction_k')
ax2.set_xlabel("Time (μs)")
ax2.set_ylabel("Forces (N)")
ax2.legend()
ax2.grid(True)

# ax2.set_ylabel("`", color='white')
# fig.text(0, 0.92, '(a)', fontsize=16, ha='left', va='top')
# # fig.text(0.535, 0.92, '(b)', fontsize=16, ha='left', va='top')
# # fig.text(0.535, 0.54, 'Forces (N)', fontsize=16, ha='left', va='center', rotation=90)
# fig.text(0.55, 0.92, '(b)', fontsize=16, ha='left', va='top')
# fig.text(0.55, 0.54, 'Forces (N)', fontsize=16, ha='left', va='center', rotation=90)
# ax2.yaxis.set_label_coords(-0.15, 0.5)

fig.suptitle(r"Engagement Time Simulation, $w_s$ = {} mm".format(int(w_pin * 1e3)), fontsize=16)

fig.text(0, 0.92, "(a)", ha='left', va='top')
fig.text(0.525, 0.92, "(b)", ha='left', va='top')
# fig.text(0.55, 0.92, "(f)", ha='left', va='top')
# fig.text(0.56, 0.92, "(d)", ha='left', va='top')

# save_data = {'data': [T, Y, dY, X, dX, events, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic]}
# np.save(save_folder + timestamp + "_wpin={:0.1f}mm".format(w_pin*1e3), save_data, allow_pickle=True)
# plt.savefig(save_folder + timestamp + "_wpin={:0.1f}mm" + ".png", dpi=300)
# # # plt.savefig("figures/" + timestamp + "_wpin={:0.1f}mm" + ".svg")
# plt.savefig(save_folder + timestamp + "_wpin={:0.1f}mm" + ".pdf")

plt.show()
