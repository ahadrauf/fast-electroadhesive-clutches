import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit, root_scalar
from datetime import datetime
import csv

now = datetime.now()
name_clarifier = "_release_time_sim_paramsweep"
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

file_loc = "../data/"
file_name = "20250103_22_13_45_release_time_sim_paramsweep_v2"
data = np.load(file_loc + file_name + ".npy", allow_pickle=True)

all_V = [150, 200, 250, 300]
w_pin_range_all, t_release_sim_w_pin, t_dielectric_range_all, t_release_sim_t_dielectric, L_range_all, t_release_sim_L, tau_1090_range_all, t_release_sim_fall_time = data

fig, axs = plt.subplots(2, 2, layout='constrained', figsize=(6, 4.5))
ax1, ax2, ax3, ax4 = axs.flat

for V, w_pin_range, t_release in zip(all_V, w_pin_range_all, t_release_sim_w_pin):
    ax1.plot(w_pin_range, t_release, label="{:d} V".format(int(V)))
for V, t_dielectric_range, t_release in zip(all_V, t_dielectric_range_all, t_release_sim_t_dielectric):
    ax2.semilogx(np.array(t_dielectric_range) * 1e6, t_release, label="{:d} V".format(int(V)))
for V, L_range, t_release in zip(all_V, L_range_all, t_release_sim_L):
    ax3.semilogx(np.array(L_range) * 1e3, t_release, label="{:d} V".format(int(V)))
for V, tau_1090_range, t_release in zip(all_V, tau_1090_range_all, t_release_sim_fall_time):
    ax4.semilogx(np.array(tau_1090_range) * 1e6, t_release, label="{:d} V".format(int(V)))

ax1.set_xlabel(r"Substrate Width (mm)")
ax1.axvline(2, c='k', ls='--')
ax1.axvline(2.5, c='k', ls='--')
ax1.axvline(3, c='k', ls='--')
ax1.axvline(4, c='k', ls='--')
ax1.axvline(5, c='k', ls='--')
ax1.axvline(6, c='k', ls='--')
ax1.set_yticks([0, 5, 10])
ax1.grid(True)

ax2.set_xlabel(r"Dielectric Thickness (μm)")
ax2.axvline(24e-6 * 1e6, c='k', ls='--')
ax2.set_xticks([1, 10, 1e2, 1e3, 1e4], [r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$'])
ax2.set_yticks([0.8, 1, 1.2, 1.4])
ax2.grid(True)

ax3.set_xlabel(r"Substrate Length (mm)")
ax3.axvline(55.5e-3 * 1e3, c='k', ls='--')
ax3.set_xticks([1, 10, 1e2, 1e3], [r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$'])
ax3.grid(True)

ax4.set_xlabel(r"Fall Time $t_{90\%-10\%}$ (µs)")
ax4.set_xticks([1, 10, 1e2, 1e3, 1e4], [r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$'])
ax4.grid(True)
ax4.axvline(5.2, c='k', ls='--')
ax2.legend(fontsize=13)

fig.text(0.04, .99, "(b.i)", transform=fig.transFigure, horizontalalignment='left', verticalalignment='top')
fig.text(0.51, .99, "(b.ii)", transform=fig.transFigure, horizontalalignment='left', verticalalignment='top')
fig.text(0.04, .5, "(b.iii)", transform=fig.transFigure, horizontalalignment='left', verticalalignment='top')
fig.text(0.51, .5, "(b.iv)", transform=fig.transFigure, horizontalalignment='left', verticalalignment='top')

fig.supylabel("Release Time (ms)", fontsize=16, y=0.55)

# plt.savefig(save_folder + timestamp + ".png", dpi=300)
# # # plt.savefig("figures/" + timestamp + ".svg")
# plt.savefig(save_folder + timestamp + ".pdf")

plt.show()
