import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit, root_scalar
from datetime import datetime
import csv

now = datetime.now()
name_clarifier = "_engagement_time_sim_paramsweep"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
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
file_name = "20241207_21_42_39_engagement_time_sim_paramsweep_v2"
data = np.load(file_loc + file_name + ".npy", allow_pickle=True)

w_pin_range, t_engage_w_pin, t_dielectric_range, t_engage_t_dielectric, L_range, t_engage_L, tau_range, t_engage_tau = data
all_V = np.array([150, 200, 250, 300])
print(w_pin_range)

fig, axs = plt.subplots(2, 2, layout='constrained', figsize=(6, 4.5))  # figsize=(10, 9))
ax1, ax2, ax3, ax4 = axs.flat

w_pin_range = np.array(w_pin_range)
L_range = np.array(L_range)
for V, t_engage in zip(all_V, t_engage_w_pin):
    ax1.plot(w_pin_range * 1e3, t_engage, label="{:d} V".format(int(V)))
for V, t_engage in zip(all_V, t_engage_t_dielectric):
    ax2.semilogx(t_dielectric_range * 1e6, t_engage, label="{:d} V".format(int(V)))
for V, t_engage in zip(all_V, t_engage_L):
    ax3.semilogx(L_range * 1e3, t_engage, label="{:d} V".format(int(V)))
for V, t_engage in zip(all_V, t_engage_tau):
    ax4.semilogx(tau_range * 1e6, t_engage, label="{:d} V".format(int(V)))

ax1.set_xlabel(r"Substrate Width (mm)")
ax1.axvline(2, c='k', ls='--')
ax1.axvline(2.5, c='k', ls='--')
ax1.axvline(3, c='k', ls='--')
ax1.axvline(4, c='k', ls='--')
ax1.axvline(5, c='k', ls='--')
ax1.axvline(6, c='k', ls='--')
ax1.set_yticks([0, 300, 600])
ax1.grid(True)

ax2.set_xlabel(r"Dielectric Thickness (μm)")
ax2.axvline(24e-6 * 1e6, c='k', ls='--')
ax2.set_xticks([1, 10, 1e2, 1e3, 1e4], [r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$'])
ax2.set_yticks([0, 5, 10, 15])
ax2.grid(True)

ax3.set_xlabel(r"Substrate Length (mm)")
ax3.axvline(55.5e-3 * 1e3, c='k', ls='--')
ax3.set_xticks([1, 10, 1e2, 1e3], [r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$'])
ax3.set_yticks([0, 2, 4])
ax3.grid(True)

ax4.set_xlabel(r"Rise Time $t_{10\%-90\%}$ (µs)")
ax4.set_xticks([1, 10, 1e2, 1e3, 1e4], [r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$'])
ax4.grid(True)
ax4.set_yticks([0, 150, 300])
ax4.axvline(8.2, c='k', ls='--')
ax4.legend(fontsize=13)

fig.text(0.04, .99, "(a.i)", transform=fig.transFigure, horizontalalignment='left', verticalalignment='top')
fig.text(0.51, .99, "(a.ii)", transform=fig.transFigure, horizontalalignment='left', verticalalignment='top')
fig.text(0.04, .5, "(a.iii)", transform=fig.transFigure, horizontalalignment='left', verticalalignment='top')
fig.text(0.51, .5, "(a.iv)", transform=fig.transFigure, horizontalalignment='left', verticalalignment='top')

fig.supylabel("Engagement Time (µs)", fontsize=16, y=0.55)

# plt.savefig(save_folder + timestamp + ".png", dpi=300)
# # # plt.savefig("figures/" + timestamp + ".svg")
# plt.savefig(save_folder + timestamp + ".pdf")

plt.show()
