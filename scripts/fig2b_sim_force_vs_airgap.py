import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

now = datetime.now()
name_clarifier = "_sim_force_vs_airgap"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
plt.rcParams["font.family"] = "Arial"
plt.rc('font', size=16)  # controls default text size
plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=16)  # fontsize of the y tick labels
plt.rc('legend', fontsize=16)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
save_folder = "../figures/"

N = 4
eps0 = 8.854e-12
k = 50
k3 = 3
w_trace = 1.5e-3
s_trace = 0.5e-3
depth = 2e-3
t_dielectric = 24e-6
V = 300

t_air = np.array([1e-08, 5e-08, 1e-07, 5e-07, 1e-06, 1.5e-06, 2e-06, 2.5e-06, 3e-06, 3.5e-06, 4e-06, 4.5e-06, 5e-06, 5.5e-06, 6e-06, 6.5e-06, 7e-06, 7.5e-06, 8e-06, 8.5e-06, 9e-06, 9.5e-06, 1e-05])
F_sim = np.array([4.983213165, 4.259480874, 3.557678567, 1.251076374, 0.551608166, 0.30944998,
                  0.198046354, 0.137684524, 0.10129898, 0.077684681, 0.061488456, 0.049895612,
                  0.041310952, 0.034770727, 0.029674037, 0.025629321, 0.022360677, 0.019683422,
                  0.017462908, 0.015598693, 0.014019571, 0.012669457, 0.011505501])
t_air_k3 = np.array([1e-08, 1e-7, 5e-7, 1e-6, 1.5e-6, 2e-6, 2.5e-06, 3e-06, 3.5e-06,
                     4e-06, 4.5e-06, 5e-06, 5.5e-06,
                     6e-06, 6.5e-06, 7e-06, 7.5e-06,
                     8e-06, 8.5e-06, 9e-06, 9.5e-06,
                     1e-05])
F_sim_k3 = np.array([0.018698892813421, 0.018288348708032, 0.016609949321657, 0.014819386358697, 0.013304825292561, 0.012012269988469, 0.010899977174276, 0.009936122400096, 0.009093841164941,
                     0.00835576803779, 0.007703809700956, 0.007126348391546, 0.006610616718979,
                     0.006149831043665, 0.0057351863742, 0.005360688717341, 0.005023272622674,
                     0.004716873701872, 0.004437426287405, 0.00418155544273, 0.003947764445824,
                     0.003732743279137])

F_sim /= depth * (N * (w_trace + s_trace)) * 1e3  # N/mm^2 --> kPa
t_air_est = np.linspace(1e-8, 1e-5, 1000)
# F_est = 1 / 8 * eps0 * (k**2) * (V**2) * (N * (w_trace + t_dielectric) * (depth + t_dielectric)) / (t_dielectric + k * t_air_est)**2  # fringing field correction - not included in paper
F_est = 1 / 8 * eps0 * (k**2) * (V**2) * (N * w_trace * depth) / (t_dielectric + k * t_air_est)**2
F_est /= depth * (N * (w_trace + s_trace)) * 1e3
# F_est_noairgap = 1 / 8 * eps0 * k * (V**2) * (N * (w_trace + t_dielectric) * (depth + t_dielectric)) / t_dielectric**2 * np.ones_like(t_air_est)  # fringing field correction - not included in paper
F_est_noairgap = 1 / 8 * eps0 * k * (V**2) * (N * w_trace * depth) / t_dielectric**2 * np.ones_like(t_air_est)
F_est_noairgap /= depth * (N * (w_trace + s_trace)) * 1e3

F_sim_k3 /= depth * (N * (w_trace + s_trace)) * 1e3  # N/mm^2 --> kPa
F_est_k3 = 1 / 8 * eps0 * (k3**2) * (V**2) * (N * w_trace * depth) / (t_dielectric + k3 * t_air_est)**2
F_est_k3 /= depth * (N * (w_trace + s_trace)) * 1e3
F_est_noairgap_k3 = 1 / 8 * eps0 * k3 * (V**2) * (N * w_trace * depth) / t_dielectric**2 * np.ones_like(t_air_est)
F_est_noairgap_k3 /= depth * (N * (w_trace + s_trace)) * 1e3

colors = ["#1964B0", "#7BB0DF", "#4DB264", "#DB060B"]
fig, ax = plt.subplots(1, 1, layout='constrained', figsize=(6, 4))
data1 = plt.scatter(t_air * 1e6, F_sim, zorder=99, c='k')
est1, = plt.plot(t_air_est * 1e6, F_est, c=colors[0], lw=2)
estdea1, = plt.plot(t_air_est * 1e6, F_est_noairgap, c=colors[0], lw=2, ls='--')

data2 = plt.scatter(t_air_k3 * 1e6, F_sim_k3, zorder=99, c='k', marker='s')
est2, = plt.plot(t_air_est * 1e6, F_est_k3, c=colors[3], lw=2)
estdea2, = plt.plot(t_air_est * 1e6, F_est_noairgap_k3, c=colors[3], lw=2, ls='--')

plt.grid(True)
plt.xlabel("Air Gap (µm)")
plt.ylabel("Normal Pressure (kPa)")
plt.legend([r"Simulated ($\kappa = 50$)", r"Model with Air Gap (Eq. 4, $\kappa = 50$)", r"Model with Zero Air Gap (Eq. 5, $\kappa = 50$)",
            r"Simulated ($\kappa = 3$)", r"Model with Air Gap (Eq. 4, $\kappa = 3$)", r"Model with Zero Air Gap (Eq. 5, $\kappa = 3$)"],
           fontsize=11.5)
plt.yscale("log")

# plt.savefig(save_folder + timestamp + ".png", dpi=300)
# plt.savefig(save_folder + timestamp + ".svg")
# plt.savefig(save_folder + timestamp + ".pdf")
plt.show()
