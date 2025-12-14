import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit, root_scalar, minimize, NonlinearConstraint
from scipy.integrate import quad
from scipy.special import kv
from datetime import datetime
from sim_release_time import *

now = datetime.now()
name_clarifier = "_force_capacity_vs_voltage_frequency_newwheatstonedata"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
save_folder = "../data/"

plt.rcParams["font.family"] = "Arial"
plt.rc('font', size=16)  # controls default text size
plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=16)  # fontsize of the y tick labels
plt.rc('legend', fontsize=16)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

markers = ['o', '^', 'D', 'x', 's', '2', '*', '3']
labels = ["DC", r"$10^0$ Hz", r"$10^1$ Hz", r"$10^2$ Hz", r"$10^3$ Hz", r"$2 \cdot 10^3$ Hz", r"$5 \cdot 10^3$ Hz", r"$10^4$ Hz", r"$10^4$ Hz (Model)"]
colors = list(reversed(["#882D71", "#1964B0", "#7BB0DF", "#4DB264", "#EE8026", "#90C987", "#71190D", "#DB060B"]))  # "#F7F057" = yellow

####################################################################################################
# Calculate Shear Force
N = 28
eps0 = 8.854e-12
w_trace = 1.5e-3
s_trace = 0.5e-3
L_pin = 95e-3
g = 9.81  # m/s^2
w_pin = 2e-3
depth_pin = 2e-3
# L_dielectric = 55.5e-3
t_dielectric = 24e-6
# w_dielectric = w_pin

V = 250
frequency = 1000
F_preload = 0.125
Fes_scalar = 1.70824E-05 * (V**2) - 0.01352372 * V + 3.208127237  # 1.66346E-05 * (V**2) - 0.013299928 * V + 3.163470918  # -0.005814336 * V + 2.373325164

mu_objet_s = 0.281  # 0.173  # 0.28
mu_pvdf_s = 0.188  # 0.154  # 0.38
mu_objet_k = 0.173
mu_pvdf_k = 0.154
mu_tot_k = mu_objet_k + mu_pvdf_k
mu_tot_s = mu_objet_s + mu_pvdf_s
# e_cc_max = 54.217
e_high, e_low, tau, alpha = 4, 54.2, 2.82e-6, 0.438
e_cc_max = np.real(e_high + (e_low - e_high) / (1 + np.power(1j * 2 * np.pi * frequency * tau, 1 - alpha)))

h2over4 = lambda h: np.square(h) / 4
h2 = lambda h: np.square(h)
F32 = lambda h: np.exp(-h2over4(h)) * np.sqrt(h) * ((h2(h) + 1) * kv(0.25, h2over4(h)) - h2(h) * kv(0.75, h2over4(h))) / 4 / np.sqrt(np.pi)

sigma_gw, C_gw = 2.76314774e-06, 3.05500533e+14
F_es = lambda k, t_air, V: Fes_scalar / 8 * eps0 * np.square(k) * (V**2) * (N * w_trace * w_pin) / np.sqrt(2 * np.pi) / sigma_gw * \
                           quad(lambda y: np.divide(np.exp(-np.square(t_air - y) / 2 / np.square(sigma_gw)), np.square(t_dielectric + k * y)), -t_dielectric / k * 0.83, 1e-3)[0]
F_gw = lambda h: C_gw * np.power(sigma_gw, 1.5) * (N * w_trace * w_pin) * F32(h / sigma_gw)  # L_dielectric


####################################################################################################

def calculate_shear_force_pred(w_dielectric, L_dielectric, F_preload):
    m_pvdf = 1.78e3 * t_dielectric * L_dielectric * w_dielectric

    t_air_with_Fes = root_scalar(lambda t_air: F_gw(t_air) - F_es(k=e_cc_max, t_air=t_air, V=V) - m_pvdf * g - F_preload, bracket=(1e-6, 10e-6)).root
    t_air_without_Fes = root_scalar(lambda t_air: F_gw(t_air) - m_pvdf * g - F_preload, x0=3e-6, x1=4e-6).root

    return (mu_pvdf_s * F_gw(t_air_with_Fes) - mu_pvdf_k * F_gw(t_air_without_Fes)) / (w_pin * L_dielectric) / 1e3


# w_dielectric = 2
# L_dielectric = 55.5
# Fshear = calculate_shear_force_pred(w_dielectric=w_dielectric, L_dielectric=L_dielectric)
# _, _, _, _, _, t_release, t_release_10pct, _, _, _, _, _, _ = sim_release_time(L_dielectric * 1e-3, w_dielectric * 1e-3, depth_pin, V, t_dielectric, frequency, F_preload, loadcell_start_ratio=0.8)
# print(Fshear, t_release)

Fshear_min = 5
# x[0] = w_dielectric, x[1] = L_dielectric, x[2] = F_preload
t_release = lambda x: sim_release_time(x[1] * 1e-3, x[0] * 1e-3, depth_pin, V, t_dielectric, frequency, F_preload, loadcell_start_ratio=0.8)[5]
# t_release = lambda x: sim_release_time(x[1] * 1e-3, x[0] * 1e-3, depth_pin, V, t_dielectric, frequency, x[2], loadcell_start_ratio=0.8)[5]
Fshear = lambda x: calculate_shear_force_pred(w_dielectric=x[0] * 1e-3, L_dielectric=x[1] * 1e-3, F_preload=F_preload) - Fshear_min
# x0 = np.array([2, 55.5])
# print(t_release(x0), Fshear(x0))

x0 = np.array([2., 30.])

for method in ["COBYLA"]:  # , "SLSQP", "trust-constr"]:
    start_time = datetime.now()

    print("Starting method:", method)
    # sol = minimize(fun=t_release, x0=x0, method=method, bounds=((1, 20), (1, 100), (0.001, 2)),
    sol = minimize(fun=t_release, x0=x0, method=method, bounds=((1., 20.), (1., 100.)),
                   constraints={'type': 'ineq', 'fun': Fshear})
    # bounds = Nelder-Mead, L-BFGS-B, TNC, SLSQP, Powell, trust-constr, COBYLA, and COBYQA
    # constraints = COBYLA, COBYQA, SLSQP and trust-constr
    print(sol)
    print("Fshear", Fshear(sol.x))
    end_time = datetime.now()
    print("Time Elapsed:", (end_time - start_time).total_seconds(), "seconds")

# fig2, axs = plt.subplots(1, 2, figsize=(12, 4.5), layout='constrained')
# ax1, ax2 = axs
# ax1_lines = []
#
# # np.save(save_folder + timestamp + "_plotdata",  [all_V, all_freq, all_x1, all_y1, all_yerr1, all_x2, all_y2, all_yerr2,
# #                                                  all_load_cell_disengage_forces, all_load_cell_preload_forces])
# # plt.savefig(save_folder + timestamp + ".png", dpi=300)
# # # # plt.savefig("figures/" + timestamp + ".svg")
# # plt.savefig(save_folder + timestamp + ".pdf")
#
# plt.show()
