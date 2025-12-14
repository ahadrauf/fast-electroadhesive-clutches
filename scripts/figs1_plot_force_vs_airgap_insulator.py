import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd

now = datetime.now()
name_clarifier = "_sim_force_vs_kappa"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
plt.rcParams["font.family"] = "Arial"
plt.rc('font', size=16)  # controls default text size
plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=16)  # fontsize of the y tick labels
plt.rc('legend', fontsize=16)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
save_folder = "C:/Users/ahadrauf/Desktop/Research/latex/electroadhesive_dynamics_paper/figures/"

N = 4
eps0 = 8.854e-12
k = 50
w_trace = 1.5e-3
s_trace = 0.5e-3
depth = 2e-3
t_dielectric = 24e-6
V = 300

kappa_range = np.power(10, np.linspace(0, 3, 12))
kappa_range = np.array(list(kappa_range) + [kappa_range[-1] * (kappa_range[1] / kappa_range[0])])
print(kappa_range)

force_sim = np.array([1.03731362125E-9, 1.01047715918E-5, 4.22170122455E-5, 1.3988422584E-4, 4.17459166827E-4, 0.001137191394394,
                      0.002793742586786, 0.006051326335745, 0.011274964541779, 0.017897686693952, 0.02460799920186, 0.030293885430629,
                      0.041307356773119])

data = np.array(pd.read_excel('../data/20250103 - FEA Insulator Voltage on Dielectric Bottom vs. Substrate Bottom.xlsx',
                              sheet_name="vs_kappa"))

colors = ["#AE76A3", "#882E72", "#1965B0", "#5289C7", "#7BAFDE",  # "#D1BBD7" = light purple (hard to see)
          "#4EB265", "#90C987", "#CAE0AB", "#F6C141", "#F1932D",  # "#F7F056" = yellow
          "#E8601C", "#DC050C", "#777777"]
fig, axs = plt.subplots(1, 2, figsize=(12, 4), layout='constrained')
ax1, ax2 = axs.flat
for i in range(13):
    x_dielectric = data[1:len(data) // 2 + 1, 2 * i] - 4e-3
    V_dielectric = data[1:len(data) // 2 + 1, 2 * i + 1]
    x_insulator = (data[len(data) // 2, 0] + (data[len(data) // 2 + 1, 2 * i] - data[len(data) // 2 + 1:, 2 * i]))[::-1] - 4e-3
    V_insulator = data[len(data) // 2 + 1:, 2 * i + 1][::-1]
    idx = np.where((-2e-3 <= x_dielectric) & (x_dielectric <= 2e-3))
    x_dielectric = x_dielectric[idx]
    V_dielectric = V_dielectric[idx]
    x_insulator = x_insulator[idx]
    V_insulator = V_insulator[idx]
    print(i, kappa_range[i], "Max gap voltage:", np.max(V_dielectric - V_insulator), np.max(V_insulator) - np.min(V_insulator),
          "Average gap voltage:", np.mean(np.abs(V_dielectric - V_insulator)), np.mean(V_insulator) - np.mean(V_insulator))

    ls = '--' if i == 12 else '-'
    ax1.plot(x_dielectric * 1e3, V_dielectric - V_insulator, label=r"$T_{air}$ = " + "{} μm".format(i + 1),
             color=colors[i], ls=ls)

ax1.set_xlabel("X Position (mm)")
ax1.set_ylabel(r"$V_{dielectric} - V_{substrate}$ (V)")
ax1.grid(True)
ax1.annotate(r"Increasing $\kappa_s$", xy=(-1, 160), xycoords='data', xytext=(-1, -30),
             textcoords='data', va='center', ha='center', fontsize=16,
             arrowprops=dict(facecolor='k', arrowstyle='-|>', connectionstyle='arc3', relpos=(0.5, 0.5), lw=2),
             bbox=dict(pad=3, visible=False, fill=False, fc=None, ec=None),
             color='k')

force_sim /= depth * (N * (w_trace + s_trace)) * 1e3  # N/mm^2 --> kPa
ax2.plot(kappa_range, force_sim, color='k', marker='o')
ax2.grid(True)
ax2.set_xscale('log')
ax2.set_xlabel(r"Substrate Dielectric Constant $\kappa_s$")
ax2.set_ylabel("Normal Pressure (kPa)")
ax2.set_xticks([1, 10, 100, 1000, kappa_range[-1]], [r"$10^0$", r"$10^1$", r"$10^2$", r"$10^3$", r"$\infty$"])
fig.text(0, 1, "(c)", transform=fig.transFigure, horizontalalignment='left', verticalalignment='top')
fig.text(0.51, 1, "(d)", transform=fig.transFigure, horizontalalignment='left', verticalalignment='top')

# plt.savefig("../TMech 2024/figures_test/" + timestamp + ".png", dpi=300)
# # plt.savefig(save_folder + timestamp + ".svg")
# plt.savefig("../TMech 2024/figures_test/" + timestamp + ".pdf")
plt.show()

