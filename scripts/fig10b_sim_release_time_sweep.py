import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit, root_scalar
from scipy.integrate import solve_ivp
from scipy.special import kv
from datetime import datetime
from sim_release_time import *

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

file_loc = "../data/strain_tests/"

# Plot data
colors = ["#1964B0", "#F1932D", "#4DB264", "#DB060B"]  # "#F7F057" = yellow
fig, axs = plt.subplots(2, 2, layout='constrained', figsize=(6.4 * 1.3, 4.8 * 1.3))  # figsize=(10, 9))
ax1, ax2 = axs[0]
ax3, ax4 = axs[1]
L_dielectric_nom, w_pin_nom, depth_pin_nom, k_dielectric_nom, t_dielectric_nom, Fext_nom, f_nom = 55.5e-3, 2e-3, 2e-3, 54.2, 24e-6, 0.125, 1000  # F_extnom = 0.38
loadcell_start_ratio = 0.8
period_nom = 1 / 2 / f_nom
V_range = np.linspace(150, 300, 4)


def poly_fit(x, a, n, b):
    return a * np.power(x, n) + b


num_sim = 50
w_pin_range = np.power(10, np.linspace(-3, -1, num_sim))
t_release_w_pin = []
for V in V_range:
    t_release = []
    for w_pin in w_pin_range:
        # print(L_dielectric_nom, w_pin, depth_pin_nom, V, t_dielectric_nom, f_nom, Fext_nom, loadcell_start_ratio)
        T, Y, dY, X, dX, events, t_release_10pct, _, _, _, _, _, _ = sim_release_time(L_dielectric=L_dielectric_nom,
                                                                                      w_pin=w_pin, depth_pin=depth_pin_nom,
                                                                                      V_max=V, t_dielectric=t_dielectric_nom,
                                                                                      frequency=f_nom, Fext=Fext_nom,
                                                                                      loadcell_start_ratio=loadcell_start_ratio)
        t_release.append(events)
        print((datetime.now() - start_time).total_seconds() / 60, "minutes have elapsed, w_pin = {}, t_release = {}".format(w_pin, events))

    ax1.plot(w_pin_range * 1e3, t_release, label="{:d} V".format(int(V)))
    t_release_w_pin.append(t_release)
    print("Results for V =", V)
    print(w_pin_range)
    print(t_release)

ax1.set_xlabel("Substrate Width (mm)")
ax1.axvline(2, c='k', ls='--')
ax1.axvline(2.5, c='k', ls='--')
ax1.axvline(3, c='k', ls='--')
ax1.axvline(4, c='k', ls='--')
ax1.axvline(5, c='k', ls='--')
ax1.axvline(6, c='k', ls='--')
ax1.set_xticks([1, 50, 100])
ax1.grid(True)
print("Results for w_pin:")
print(w_pin_range)
print(t_release_w_pin)
print("Total runtime:", (datetime.now() - start_time).total_seconds(), 'seconds, or', (datetime.now() - start_time).total_seconds() / 60, "minutes")

t_dielectric_range = np.power(10, np.linspace(-6, -2, num_sim))
t_dielectric_range_all = []
t_release_t_dielectric = []
for V in V_range:
    t_release = []
    t_dielectric_range_curr = []
    for t_dielectric in t_dielectric_range:
        try:
            T, Y, dY, X, dX, events, t_release_10pct, _, _, _, _, _, _ = sim_release_time(L_dielectric=L_dielectric_nom,
                                                                                          w_pin=w_pin_nom, depth_pin=depth_pin_nom,
                                                                                          V_max=V, t_dielectric=t_dielectric,
                                                                                          frequency=f_nom, Fext=Fext_nom,
                                                                                          loadcell_start_ratio=loadcell_start_ratio)
            t_dielectric_range_curr.append(t_dielectric)
            t_release.append(events)
        except Exception as e:
            print("Failed sim, continuing on")

    ax2.semilogx(np.array(t_dielectric_range_curr) * 1e6, t_release)
    t_release_t_dielectric.append(t_release)
    t_dielectric_range_all.append(t_dielectric_range_curr)
    print("Results for t_dielectric, V =", V)
    print(t_dielectric_range_curr)
    print(t_release)
ax2.set_xlabel("Dielectric Thickness (μm)")
ax2.axvline(t_dielectric_nom * 1e6, c='k', ls='--')
ax2.set_xticks([1, 10, 1e2, 1e3, 1e4], [r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$'])
ax2.tick_params(axis='x', which='minor', length=2, width=1, bottom=True)
ax2.grid(True)
ax2.legend(fontsize=16)
print("Results for t_dielectric:")
print(t_dielectric_range)
print(t_release_t_dielectric)
print("Total runtime:", (datetime.now() - start_time).total_seconds(), 'seconds, or', (datetime.now() - start_time).total_seconds() / 60, "minutes")

L_range = np.power(10, np.linspace(-3, 0, num_sim))
t_release_L = []
for V in V_range:
    t_release = []
    for L in L_range:
        T, Y, dY, X, dX, events, t_release_10pct, _, _, _, _, _, _ = sim_release_time(L_dielectric=L,
                                                                                      w_pin=w_pin_nom, depth_pin=depth_pin_nom,
                                                                                      V_max=V, t_dielectric=t_dielectric_nom,
                                                                                      frequency=f_nom, Fext=Fext_nom,
                                                                                      loadcell_start_ratio=loadcell_start_ratio)
        t_release.append(events)

    ax3.semilogx(L_range * 1e3, t_release, label="{:d} V".format(int(V)))
    t_release_L.append(t_release)
ax3.set_xlabel(r"Substrate Length (mm)")
ax3.axvline(L_dielectric_nom * 1e3, c='k', ls='--')
ax3.set_xticks([1, 10, 1e2, 1e3], [r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$'])
ax3.grid(True)
print("Results for L:")
print(L_range)
print(t_release_L)
print("Total runtime:", (datetime.now() - start_time).total_seconds(), 'seconds, or', (datetime.now() - start_time).total_seconds() / 60, "minutes")

tau_1090_range = np.power(10, np.linspace(-6, -2, num_sim))
t_release_Fext = []
for V in V_range:
    t_release = []
    for tau_1090 in tau_1090_range:
        T, Y, dY, X, dX, events, t_release_10pct, _, _, _, _, _, _ = sim_release_time(L_dielectric=L,
                                                                                      w_pin=w_pin_nom, depth_pin=depth_pin_nom,
                                                                                      V_max=V, t_dielectric=t_dielectric_nom,
                                                                                      frequency=f_nom, Fext=Fext_nom,
                                                                                      loadcell_start_ratio=loadcell_start_ratio,
                                                                                      tau_1090=tau_1090)
        t_release.append(events)
    ax4.semilogx(tau_1090_range * 1e6, t_release, label="{:d} V".format(int(V)))
    t_release_Fext.append(t_release)
ax4.set_xlabel(r"Fall Time $t_{90\%-10\%}$ (µs)")
ax4.set_xticks([1, 10, 1e2, 1e3, 1e4], [r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$'])
ax4.grid(True)

fig.supylabel("Release Time (ms)", fontsize=20)

print(timestamp)
# np.save(save_folder + timestamp, [w_pin_range, t_release_w_pin, t_dielectric_range, t_release_t_dielectric, L_range, t_release_L, tau_1090_range, t_release_Fext], allow_pickle=True)

# plt.savefig(save_folder + timestamp + ".png", dpi=300)
# # # plt.savefig("figures/" + timestamp + ".svg")
# plt.savefig(save_folder + timestamp + ".pdf")

end_time = datetime.now()
print("Total runtime:", (end_time - start_time).total_seconds(), 'seconds, or', (end_time - start_time).total_seconds() / 60, "minutes")
plt.show()
