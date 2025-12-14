import numpy as np
import matplotlib.pyplot as plt
from pandas import read_excel
from datetime import datetime
from scipy.optimize import curve_fit
from scipy.special import kv
import re

now = datetime.now()
name_clarifier = "_lcr_capacitance_vs_force"
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

folder = "../data/"
file_name = "20240929 - IDF Capacitance vs. Load Force Vertical Test Setup.xlsx"

df_sheet_index = read_excel(folder + file_name, sheet_name=None)
fig, axs = plt.subplots(1, 2, layout='constrained', figsize=(6*4/2.75, 4))
ax1, ax2 = axs

remove_nan = lambda x: x[~np.isnan(x)]
float_template = r"([+-]?([0-9]+([.][0-9]*)?|[.][0-9]+))"
forces, caps, gaps = {}, {}, {}
for sheet_name in df_sheet_index.keys():
    data = df_sheet_index[sheet_name].to_numpy()
    match = re.search(r"w=" + float_template, sheet_name)
    width = float(match.group(1))
    if width == 6.35:
        continue
    if width not in forces:
        forces[width] = []
        caps[width] = []
        gaps[width] = []
    forces[width].append(remove_nan(np.array(data[:, 0], dtype='float')))
    caps[width].append(remove_nan(np.array(data[:, 1], dtype='float')))
    gaps[width].append(remove_nan(np.array(data[:, 2], dtype='float')))

colors = ["#882E72", "#1965B0", "#4EB265", "#E8601C"]
h2over4 = lambda h: np.square(h)/4
h2 = lambda h: np.square(h)
F32 = lambda h: np.exp(-h2over4(h))*np.sqrt(h)/4/np.sqrt(np.pi)*((h2(h) + 1)*kv(0.25, h2over4(h)) - h2(h)*kv(0.75, h2over4(h)))
F52 = lambda h: np.exp(-h2over4(h))*np.power(h, 1.5)*((2*h**2 + 3)*kv(0.75, h2over4(h)) - (2*np.square(h) + 5)*kv(0.25, h2over4(h)))/8/np.sqrt(np.pi)

length = 60e-3
avg_fit_albion = []
for idx, width in enumerate(forces.keys()):
    print("Width =", width)
    f_range = np.linspace(1e-6, 8, 1000)
    caps_interp = np.zeros((len(f_range), len(forces[width])))
    gaps_interp = np.zeros((len(f_range), len(forces[width])))
    for i in range(len(forces[width])):
        caps_interp[:, i] = np.interp(f_range, forces[width][i], caps[width][i])
        gaps_interp[:, i] = np.interp(f_range, forces[width][i], gaps[width][i])
    avg_cap, std_cap = np.nanmean(caps_interp, axis=1), np.nanstd(caps_interp, axis=1)
    avg_gap, std_gap = np.nanmean(gaps_interp, axis=1), np.nanstd(gaps_interp, axis=1)
    ax1.plot(f_range, avg_cap, label=r"$w_s =$" + "{} mm".format(width), c=colors[idx])
    ax1.fill_between(f_range, avg_cap - std_cap, avg_cap + std_cap, alpha=0.3, color=colors[idx])
    ax2.plot(f_range, avg_gap, label=r"$w_s =$" + "{} mm".format(width), c=colors[idx])
    ax2.fill_between(f_range, avg_gap - std_gap, avg_gap + std_gap, alpha=0.3, color=colors[idx])

    # Measure contact stiffness
    idx_7N = np.argmin(np.abs(f_range - 5))
    sol = np.polyfit(avg_gap[idx_7N:], f_range[idx_7N:], 1)
    print("Stiffness:", sol)

    F_gw = lambda h, sigma, C: C*np.power(sigma, 1.5)*(width*1e-3)*length*F32(h/sigma)
    # F_gt = lambda h, sigma, C: C * np.power(sigma, 2.5) * (width * 1e-3) * length * F52(h / sigma)
    avg_gap_fit = avg_gap[:len(avg_gap)]
    f_range_fit = f_range[:len(avg_gap)]
    print("Gap to fit from", avg_gap[len(avg_gap)//4])
    sol = curve_fit(F_gw, avg_gap_fit/1e6, f_range_fit, p0=(5e-6, 2e15), maxfev=10000)[0]  # , bounds=((0, 0), (1e-5, 10000)))[0]

    gap_sim = [avg_gap[0]]
    dgap_sim = (avg_gap[0] - avg_gap[-1])/1000
    F_pred = [F_gw(avg_gap[0]/1e6, *sol)]
    while F_pred[-1] < 8:
        curr_gap = avg_gap[0] - dgap_sim*len(gap_sim)
        gap_sim.append(curr_gap)
        F_pred.append(F_gw(curr_gap/1e6, *sol))

    if 2 <= width <= 3:
        avg_fit_albion.append((sol))
    C_sim = 28*45*8.854e-12*(width*1e-3*1.5e-3)/4/(24e-6 + 45*1e-6*np.array(gap_sim[:-1]))
    if idx == 2:
        ax1.plot(F_pred[:-1], C_sim*1e12, ls='--', c=colors[idx], label="Model")
        ax2.plot(F_pred[:-1], gap_sim[:-1], ls='--', c=colors[idx], label="Model")
    else:
        ax1.plot(F_pred[:-1], C_sim*1e12, ls='--', c=colors[idx])
        ax2.plot(F_pred[:-1], gap_sim[:-1], ls='--', c=colors[idx])
    print("C_gw -->", sol)
    F_pred = F_gw(avg_gap/1e6, *sol)
    err = np.abs(np.divide(F_pred - f_range, f_range)*100)
    print("Error", np.mean(err))

Cgw_albion = np.mean(avg_fit_albion, axis=0)  # N/m^3.5
print("Cgw_albion", Cgw_albion)
print("Albion Alloys Avg. sigma_gw", Cgw_albion[0])
print("Albion Alloys Cgw: {} N/m^3.5, {} N/mm^3.5, {} mN/um^3.5".format(Cgw_albion[1],
                                                                        Cgw_albion[1]/(1000**3.5),
                                                                        Cgw_albion[1]*1000/(1000000**3.5)),
      avg_fit_albion)
print("Albion Alloys Cgw - 1 std:", Cgw_albion[1] - np.std(avg_fit_albion, axis=0)[1])

ax1.set_ylabel("Measured Capacitance (pF)")
ax2.set_ylabel("Air Gap " + r"$T_{air}$" + " (µm)")
ax1.set_xlabel("Normal Force (N)")
ax2.set_xlabel("Normal Force (N)")
ax1.set_xticks([0, 2, 4, 6, 8])
ax2.set_xticks([0, 2, 4, 6, 8])
ax1.set_yticks([10, 20, 30, 40, 50])
ax1.grid(True)
ax2.grid(True)
ax1.legend(fontsize=16)
ax2.legend(fontsize=16)

fig.text(0.0, 1, "(b)", transform=fig.transFigure, horizontalalignment='left', verticalalignment='top')
fig.text(0.5, 1, "(c)", transform=fig.transFigure, horizontalalignment='left', verticalalignment='top')

# plt.savefig(save_folder + timestamp + ".png", dpi=300)
# # # plt.savefig("figures/" + timestamp + ".svg")
# plt.savefig(save_folder + timestamp + ".pdf")
plt.show()
