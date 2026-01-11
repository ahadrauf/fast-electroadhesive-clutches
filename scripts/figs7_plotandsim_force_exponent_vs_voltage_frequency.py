import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit, root_scalar
from scipy.integrate import quad
from scipy.special import kv
from datetime import datetime

now = datetime.now()
name_clarifier = "_force_capacity_vs_voltage_frequency_newwheatstonedata"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
save_folder = "../data/"

plt.rcParams["font.family"] = "Arial"
plt.rc('font', size=24)  # controls default text size
plt.rc('axes', labelsize=20)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=20)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=20)  # fontsize of the y tick labels
plt.rc('legend', fontsize=20)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

file_name = "20241114_14_14_33_force_capacity_vs_voltage_frequency_newwheatstonedata_plotdata"
data = np.load(save_folder + file_name + ".npy", allow_pickle=True)
all_V, all_freq, all_x1, all_y1, all_yerr1, all_x2, all_y2, all_yerr2, all_load_cell_disengage_forces, all_load_cell_preload_forces, all_load_cell_max_force = data

fig2, axs = plt.subplots(1, 2, figsize=(12, 4.5), layout='constrained')
ax1, ax2 = axs
markers = ['o', '^', 'D', 'x', 's', '2', '*', '3']
labels = ["DC", r"$10^0$ Hz", r"$10^1$ Hz", r"$10^2$ Hz", r"$10^3$ Hz", r"$2 \cdot 10^3$ Hz", r"$5 \cdot 10^3$ Hz", r"$10^4$ Hz", r"$10^4$ Hz (Model)"]
colors = list(reversed(["#882D71", "#1964B0", "#7BB0DF", "#4DB264", "#EE8026", "#90C987", "#71190D", "#DB060B"]))  # "#F7F057" = yellow
ax1_lines = []


def exponential_func(v, a, n):
    return a * np.power(v, n)


def fit_to_curve(v, f):
    if v[0] != 0:
        v = np.insert(v, 0, 0)
        f = np.insert(f, 0, 0)
    popt, pcov = curve_fit(exponential_func, v, f, bounds=([0, 0.5], [100, 4]))
    perr = np.sqrt(np.diag(pcov))
    ss_res = np.sum(np.square(f - exponential_func(v, *popt)))
    ss_tot = np.sum(np.square(f - np.mean(f)))
    r2 = 1 - ss_res / ss_tot
    # mape = 100/len(v) * np.sum(np.abs(np.divide(f - exponential_func(v, *popt), f)))
    return popt, perr, r2


def calculate_shear_force_pred(V_max, Fpreload, frequency, Fes_scalar=1.):
    N = 28
    eps0 = 8.854e-12
    w_trace = 1.5e-3
    s_trace = 0.5e-3
    L_pin = 95e-3
    g = 9.81  # m/s^2
    w_pin = 2e-3
    L_dielectric = 55.5e-3
    t_dielectric = 24e-6
    w_dielectric = w_pin
    # Fpreload *= 10

    mu_objet_s = 0.281  # 0.173  # 0.28
    mu_pvdf_s = 0.188  # 0.154  # 0.38
    mu_objet_k = 0.173
    mu_pvdf_k = 0.154
    mu_tot_k = mu_objet_k + mu_pvdf_k
    mu_tot_s = mu_objet_s + mu_pvdf_s
    # e_cc_max = 54.217
    e_high, e_low, tau, alpha = 4, 54.2, 2.82e-6, 0.438
    e_cc_max = np.real(e_high + (e_low - e_high) / (1 + np.power(1j * 2 * np.pi * frequency * tau, 1 - alpha)))
    m_pvdf = 1.78e3 * t_dielectric * L_dielectric * w_dielectric

    h2over4 = lambda h: np.square(h) / 4
    h2 = lambda h: np.square(h)
    F32 = lambda h: np.exp(-h2over4(h)) * np.sqrt(h) * ((h2(h) + 1) * kv(0.25, h2over4(h)) - h2(h) * kv(0.75, h2over4(h))) / 4 / np.sqrt(np.pi)

    sigma_gw, C_gw = 2.76314774e-06, 3.05500533e+14
    F_es = lambda k, t_air, V: Fes_scalar / 8 * eps0 * np.square(k) * (V**2) * (N * w_trace * w_pin) / np.sqrt(2 * np.pi) / sigma_gw * \
                               quad(lambda y: np.divide(np.exp(-np.square(t_air - y) / 2 / np.square(sigma_gw)), np.square(t_dielectric + k * y)), -t_dielectric / k * 0.83, 1e-3)[0]
    F_gw = lambda h: C_gw * np.power(sigma_gw, 1.5) * (N * w_trace * w_pin) * F32(h / sigma_gw)  # L_dielectric

    t_air_with_Fes = root_scalar(lambda t_air: F_gw(t_air) - F_es(k=e_cc_max, t_air=t_air, V=V_max) - m_pvdf * g - Fpreload, bracket=(1e-6, 10e-6)).root
    t_air_without_Fes = root_scalar(lambda t_air: F_gw(t_air) - m_pvdf * g - Fpreload, x0=3e-6, x1=4e-6).root

    # string_to_write = "V = {}, f = {}Hz, Fpreload = {:0.3f} --> t_air, with Fes={:0.3e}, without Fes={:0.3e} --> Fpred = {:0.3f}, Ppred = {:0.3f} kPa, C_fin = {:0.3e}, Qmax = {:0.3e}, tauRC = {:0.3e}"
    # print(string_to_write.format(V, frequency, Fpreload, t_air_with_Fes, t_air_without_Fes,
    #                              mu_pvdf_s * (F_gw(t_air_with_Fes) - F_gw(t_air_without_Fes)),
    #                              mu_pvdf_s * (F_gw(t_air_with_Fes) - F_gw(t_air_without_Fes)) / (w_pin * L_dielectric) / 1e3,
    #                              C(t_air_with_Fes), Ileakage * tau_leakage(t_air_with_Fes) * (1 - np.exp(-period / tau_leakage(t_air_with_Fes))),
    #                              Rleakage * C(t_air_with_Fes)))
    # return mu_pvdf_s * (F_gw(t_air_with_Fes)) / (w_pin * L_dielectric) / 1e3
    return (mu_pvdf_s * F_gw(t_air_with_Fes) - mu_pvdf_k * F_gw(t_air_without_Fes)) / (w_pin * L_dielectric) / 1e3


all_x1, all_y1, all_yerr1 = [], [], []
all_power_fits = []
all_power_fits_std = []
all_r2_scores = []
for i, f in enumerate(all_freq):
    y = [0]
    yerr = [0]
    ypred = []

    F_preload = np.mean([np.mean(all_load_cell_preload_forces[f][V]) for V in all_V])
    for V in all_V:
        print("f = {} Hz, {} V: N = {}, Max Force = {:0.3f} N, std = {:0.3f} N, Preload = {:0.3f}, raw = {}".format(f, V, len(all_load_cell_max_force[f][V]), np.mean(all_load_cell_max_force[f][V]),
                                                                                                                    np.std(all_load_cell_max_force[f][V]), np.mean(all_load_cell_preload_forces[f][V]),
                                                                                                                    all_load_cell_max_force[f][V]))
        y.append(np.mean(all_load_cell_max_force[f][V]))
        yerr.append(np.std(all_load_cell_max_force[f][V]))

    V_sim = np.linspace(0, 300, 25)
    for V in V_sim:
        Fes_scalar = 1.70824E-05 * (V**2) - 0.01352372 * V + 3.208127237  # 1.66346E-05 * (V**2) - 0.013299928 * V + 3.163470918  # -0.005814336 * V + 2.373325164
        print("V", V, "Fes_scalar", Fes_scalar)
        ypred.append(calculate_shear_force_pred(V_max=V, Fpreload=F_preload, frequency=f, Fes_scalar=Fes_scalar))  # 4086858.8994281082))
    ax1_line = ax1.errorbar([0] + all_V, y=y, yerr=yerr, fmt='-', lw=1, capsize=4, capthick=1, alpha=0.8, c=colors[i])
    ax1_points = ax1.scatter([0] + all_V, y, marker=markers[i], c=colors[i])
    ax1_lines.append((ax1_line, ax1_points))
    popt, perr, r2 = fit_to_curve([0] + all_V, y)
    print("Power fit for f = {}: popt = {}, r2 = {}".format(f, popt, r2))
    all_power_fits.append(popt[1])
    all_power_fits_std.append(perr[1])
    all_r2_scores.append(r2)

    if f == 10000:
        model_line, = ax1.plot(V_sim, ypred, ls='--', c=colors[i], label="Model")
        ax1_lines.append(model_line)

    all_x1.append([0] + all_V)
    all_y1.append(y)
    all_yerr1.append(yerr)

print("All power fits:", all_power_fits)
print(all_r2_scores)

ax2.errorbar(all_freq, y=all_power_fits, yerr=all_power_fits_std, fmt='-', lw=2, capsize=5, capthick=2, alpha=0.8, c='k')
ax2.scatter(all_freq, all_power_fits, marker='o', c='k')

ax1.grid(True)
ax1.set_xlabel("Voltage (V)")
ax1.set_ylabel("Shear Capacity (kPa)")
ax1.set_xticks([0, 50, 100, 150, 200, 250, 300])

ax2.grid(True, which='minor', c='0.9')
ax2.grid(True)
ax2.set_xlabel("Drive Frequency (Hz)")
ax2.set_ylabel(r"Shear Capacity $\propto V^n$")
ax2.set_xscale('log')
ax2.set_xticks([0.1, 1, 10, 100, 1000, 10000], ["DC", r"$10^0$", r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$"])

ax1.legend(ax1_lines, labels, loc='upper left', fontsize=13)
fig2.text(0.00, 1, "(a)", transform=fig2.transFigure, horizontalalignment='left', verticalalignment='top', fontsize=20)
fig2.text(0.505, 1, "(b)", transform=fig2.transFigure, horizontalalignment='left', verticalalignment='top', fontsize=20)

# np.save(save_folder + timestamp + "_plotdata",  [all_V, all_freq, all_x1, all_y1, all_yerr1, all_x2, all_y2, all_yerr2,
#                                                  all_load_cell_disengage_forces, all_load_cell_preload_forces])
# plt.savefig(save_folder + timestamp + ".png", dpi=300)
# # # plt.savefig("figures/" + timestamp + ".svg")
# plt.savefig(save_folder + timestamp + ".pdf")

plt.show()

