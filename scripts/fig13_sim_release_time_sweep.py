import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit, root_scalar
from scipy.integrate import solve_ivp
from scipy.special import kv
from datetime import datetime
import csv
# from sim_release_time_20241206 import *
from sim_release_time_sweep_20241209 import *

plt.rcParams["font.family"] = "Arial"
# plt.rcParams["font.family"] = "serif"
# plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rc('font', size=24)  # controls default text size
plt.rc('axes', labelsize=20)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=20)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=20)  # fontsize of the y tick labels
plt.rc('legend', fontsize=20)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
# plt.style.use('tableau-colorblind10')

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


# w_pin_range = np.array([1e-3, 2e-3, 3e-3, 4e-3, 5e-3, 6e-3, 7e-3])


def poly_fit(x, a, n, b):
    return a * np.power(x, n) + b


now = datetime.now()
# name_clarifier = "_release_time_sim_paramsweep_ws=91to100"
name_clarifier = "_release_time_sim_paramsweep_wpin=21to30_V=250and300"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
save_folder = "C:/Users/ahadrauf/Desktop/Research/latex/electroadhesive_dynamics_paper/figures_test/"
start_time = datetime.now()
num_sim = 5
w_pin_range = np.linspace(21e-3, 30e-3, num_sim)
# w_pin_range = np.linspace(60e-3, 100e-3, num_sim)
# w_pin_range = np.power(10, np.linspace(-3, -2, num_sim))

# [0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.01 ]
# t_engage_w_pin = [[0.7988561639241996, 0.9051603966756259, 1.0771356510893695, 1.2921352492611922, 1.5465758951482391, 1.8399474274205816, 2.169127049125032, 2.5192316385587565, 2.894384171843549, 3.283775161054685],
#                   [0.7973297967810035, 0.8881312231585241, 1.0258987140374725, 1.1883933514564193, 1.37390199318728, 1.581816561267912, 1.811168310288865, 2.061049886548248, 2.33004169548753, 2.612052355647749],
#                   [0.7968764041752278, 0.8716701030901913, 0.9825643878095783, 1.1137875089913576, 1.2616988875720054, 1.4258242677548438, 1.6053196634592313, 1.7994603146485522, 2.0076879613015923, 2.2294044369555808],
#                   [0.7953305974769523, 0.8590471797979109, 0.9532071028731206, 1.0644060655702974, 1.1908794957450477, 1.3315352392301043, 1.4849459114176697, 1.6504239081564023, 1.8273696119158822, 2.015289151530468]]
# for i, V in enumerate(V_range):
#     t_engage = t_engage_w_pin[i]
#     ax1.plot(w_pin_range * 1e3, t_engage, label="{:d} V".format(int(V)))

t_engage_w_pin = []
for V in V_range:
    if V == 150 or V == 200:
        t_engage = [0] * len(w_pin_range)
    else:
        t_engage = []
        for w_pin in w_pin_range:
            # print(L_dielectric_nom, w_pin, depth_pin_nom, V, t_dielectric_nom, f_nom, Fext_nom, loadcell_start_ratio)
            T, Y, dY, X, dX, events, t_release_10pct, _, _, _, _, _, _ = sim_release_time(L_dielectric=L_dielectric_nom,
                                                                                          w_pin=w_pin, depth_pin=depth_pin_nom,
                                                                                          V_max=V, t_dielectric=t_dielectric_nom,
                                                                                          frequency=f_nom, Fext=Fext_nom,
                                                                                          loadcell_start_ratio=loadcell_start_ratio)
            t_engage.append(events)
            print((datetime.now() - start_time).total_seconds() / 60, "minutes have elapsed, w_pin = {}, t_engage = {}".format(w_pin, events))

    # ax.plot(w_pin_range, t_engage, marker='o', label=r"$w_{pin}$" + " = {:0.1f} mm".format(w_pin * 1e3))
    ax1.plot(w_pin_range * 1e3, t_engage, label="{:d} V".format(int(V)))
    # ax1.semilogx(w_pin_range * 1e3, t_engage, label="{:d} V".format(int(V)))
    t_engage_w_pin.append(t_engage)
    print("Results for V =", V)
    print(w_pin_range)
    print(t_engage)

    # popt, pcov = curve_fit(poly_fit, w_pin_range * 1e3, t_engage, p0=(1, 3, 0))
    # print("Poly order:", popt[1], popt)
# ax.set_xlabel("Voltage (V)")
ax1.set_xlabel(r"$w_s$ (mm)")
# ax1.set_ylabel("Engagement Time (µs)")
ax1.axvline(2, c='k', ls='--')
ax1.axvline(2.5, c='k', ls='--')
ax1.axvline(3, c='k', ls='--')
ax1.axvline(4, c='k', ls='--')
ax1.axvline(5, c='k', ls='--')
ax1.axvline(6, c='k', ls='--')
# ax.semilogy(True)
ax1.grid(True)
# ax1.set_xticks([1, 10, 1e2], [r'$10^0$', r'$10^1$', r'$10^2$'])
# ax1.legend()
print("Results for w_pin:")
print(w_pin_range)
print(t_engage_w_pin)
print("Total runtime:", (datetime.now() - start_time).total_seconds(), 'seconds, or', (datetime.now() - start_time).total_seconds() / 60, "minutes")

# k_range = np.linspace(4, 100, num_sim)
# f_range = np.power(10, np.linspace(0, 5, 2 * num_sim))
# t_dielectric_range = np.power(10, np.linspace(np.log10(3e-6), -2, num_sim))  # 10000e-6
t_dielectric_range = np.linspace(3e-6, 10e-6, 10)  # 10000e-6
t_dielectric_range_all = []
t_engage_t_dielectric = []
# for V in V_range:
#     t_engage = []
#     # if V == 150:
#     #     t_dielectric_range_curr = np.power(10, np.linspace(np.log10(1.5e-6), -2, 25))  # 10000e-6
#     #     t_engage = 1e-3 * np.array([855.255, 843.890, 846.359, 856.449, 870.251,
#     #                                 883.596, 893.974, 901.531, 909.767, 913.489,
#     #                                 918.225, 921.967, 923.779, 926.924, 928.278,
#     #                                 931.343, 927.236, 921.483, 909.933, 891.326,
#     #                                 878.432, 881.049, 891.553, 906.541, 923.582])
#     # elif V == 200:
#     #     t_dielectric_range_curr = np.power(10, np.linspace(np.log10(1.5e-6), -2, 25))  # 10000e-6
#     #     t_engage = 1e-3 * np.array([1091.477, 876.801, 848.092, 843.662, 848.850,
#     #                                 859.236, 871.524, 883.198, 891.974, 898.977,
#     #                                 904.576, 909.671, 915.430, 918.933, 922.367,
#     #                                 924.254, 923.504, 918.977, 908.191, 890.868,
#     #                                 880.091, 882.037, 893.460, 908.447, 926.070])
#     # else:
#     # t_dielectric_range_curr = t_dielectric_range
#     if V == 150:
#         t_dielectric_range_curr = np.linspace(1e-6, 10e-6, 30)
#         t_engage = [0.9304779351350738, 0.8669051250513764, 0.8510055815035729, 0.8456023323359795, 0.8436788045867208, 0.8437999758360634, 0.8449333636579882, 0.846542311367788, 0.8486916979380994, 0.8509634576136541, 0.8533322461902748, 0.8555744671870076, 0.8581411741145579, 0.8604224029278565,
#                     0.8626752618396965, 0.8648664413654202, 0.866736724791305, 0.8688855292960221, 0.8707145048775117, 0.872481323988569, 0.8738288143165358, 0.8758036503021289, 0.8772913696800585, 0.8787184023323737, 0.8800334286586567, 0.881329058283359, 0.8823359655557489, 0.8836397038570115,
#                     0.8845081633407127, 0.8853875708046847]
#     elif V == 200:
#         t_dielectric_range_curr = np.linspace(1.5e-6, 10e-6, 30)
#         t_engage = [1.0914771174761426, 0.9260236112876631, 0.8830704615719104, 0.8648381405197026, 0.8553169442351367, 0.8500113511494338, 0.8468458261819871, 0.8449955592740347, 0.844007406992827, 0.8436049379928848, 0.8438089453515256, 0.8439013772237873, 0.8444357242805366,
#                     0.8450550080477477, 0.8459323641551867, 0.8468036173898362, 0.8477628812060729, 0.8488083822270753, 0.8498652875387126, 0.8508540443861695, 0.8519603265504714, 0.8529900913186695, 0.8540727332578275, 0.8551404620353279, 0.8562418016958876, 0.8573374393167886,
#                     0.8584208254432149, 0.8594053503162946, 0.8603968769711045, 0.8613301832829918]
#     else:
#         t_dielectric_range_curr = []
#         for t_dielectric in t_dielectric_range:
#             try:
#                 T, Y, dY, X, dX, events, t_release_10pct, _, _, _, _, _, _ = sim_release_time(L_dielectric=L_dielectric_nom,
#                                                                                               w_pin=w_pin_nom, depth_pin=depth_pin_nom,
#                                                                                               V_max=V, t_dielectric=t_dielectric,
#                                                                                               frequency=f_nom, Fext=Fext_nom,
#                                                                                               loadcell_start_ratio=loadcell_start_ratio)
#                 t_dielectric_range_curr.append(t_dielectric)
#                 t_engage.append(events)
#             except Exception as e:
#                 print("Failed sim, continuing on")
#
#     # ax.plot(w_pin_range, t_engage, marker='o', label=r"$w_{pin}$" + " = {:0.1f} mm".format(w_pin * 1e3))
#     # ax2.plot(k_range, t_engage, label="{:d} V".format(int(V)))
#     # ax2.semilogx(f_range, t_engage, label="{:d} V".format(int(V)))
#     ax2.semilogx(np.array(t_dielectric_range_curr) * 1e6, t_engage)
#     t_engage_t_dielectric.append(t_engage)
#     t_dielectric_range_all.append(t_dielectric_range_curr)
#     print("Results for t_dielectric, V =", V)
#     print(t_dielectric_range_curr)
#     print(t_engage)
# ax.set_xlabel("Voltage (V)")
# ax2.set_xlabel(r"$\kappa_s$ Dielectric Constant")
# ax2.set_xlabel("Drive Frequency (Hz)")
ax2.set_xlabel(r"$T_d$ (μm)")
# ax2.set_xticks([1, 10, 1e2, 1e3, 1e4], [r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$'])
# ax2.set_ylabel("Engagement Time (µs)")
# for f in [1, 10, 100, 1000, 2000, 5000, 10000]:
#     ax2.axvline(f, c='k', ls='--')
# # ax.semilogy(True)
# ax2.set_xticks([1, 10, 1e2, 1e3, 1e4, 1e5], [r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$'])
# ax2.xaxis.set_major_formatter(mtick.ScalarFormatter(useMathText=True))
# ax2.xaxis.set_minor_locator(plt.LogLocator(subs=np.arange(2, 10)))
ax2.tick_params(axis='x', which='minor', length=2, width=1, bottom=True)
ax2.grid(True)
# ax2.legend()
print("Results for t_dielectric:")
print(t_dielectric_range)
print(t_engage_t_dielectric)
print("Total runtime:", (datetime.now() - start_time).total_seconds(), 'seconds, or', (datetime.now() - start_time).total_seconds() / 60, "minutes")

# L_range = np.power(10, np.linspace(-2, 0, num_sim))
L_range = np.power(10, np.linspace(-3, 0, num_sim))
t_engage_L = []
# for V in V_range:
#     t_engage = []
#     for L in L_range:
#         # T, Y, dY, t_end, Fk, Fb, Fes, Fconstant = sim_engagement_time(L_dielectric=L, w_pin=w_pin_nom, depth_pin=depth_pin_nom, V_max=V, k_dielectric=k_dielectric_nom, t_dielectric=t_dielectric_nom,
#         #                                                               period=1, Fext=Fext_nom)
#         T, Y, dY, X, dX, events, t_release_10pct, _, _, _, _, _, _ = sim_release_time(L_dielectric=L,
#                                                                                       w_pin=w_pin_nom, depth_pin=depth_pin_nom,
#                                                                                       V_max=V, t_dielectric=t_dielectric_nom,
#                                                                                       frequency=f_nom, Fext=Fext_nom,
#                                                                                       loadcell_start_ratio=loadcell_start_ratio)
#         t_engage.append(events)
#
#     # ax.plot(w_pin_range, t_engage, marker='o', label=r"$w_{pin}$" + " = {:0.1f} mm".format(w_pin * 1e3))
#     ax3.semilogx(L_range * 1e3, t_engage, label="{:d} V".format(int(V)))
#     t_engage_L.append(t_engage)
# # ax.set_xlabel("Voltage (V)")
# ax3.set_xlabel(r"$L_s$ (mm)")
# # ax3.set_ylabel("Engagement Time (µs)")
# ax3.axvline(L_dielectric_nom * 1e3, c='k', ls='--')
# ax3.set_xticks([1, 10, 1e2, 1e3], [r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$'])
# # ax3.set_xticks([1e-1, 1, 10, 1e2, 1e3], [r'$10^{-1}$', r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$'])
# # ax.semilogy(True)
# ax3.grid(True)
# ax3.legend()
print("Results for L:")
print(L_range)
print(t_engage_L)
print("Total runtime:", (datetime.now() - start_time).total_seconds(), 'seconds, or', (datetime.now() - start_time).total_seconds() / 60, "minutes")

# t_dielectric_range = np.linspace(1e-6, 1000e-6, num_sim)  # 10000e-6
# t_dielectric_range = np.power(10, np.linspace(-6, -2, num_sim))  # 10000e-6
tau_1090_range = np.power(10, np.linspace(-6, -2, num_sim))
t_engage_Fext = []
# Fext_range = np.linspace(0, 1, num_sim)
# for V in V_range:
#     t_engage = []
#     for tau_1090 in tau_1090_range:
#         # T, Y, dY, t_end, Fk, Fb, Fes, Fconstant = sim_engagement_time(L_dielectric=L_dielectric_nom, w_pin=w_pin_nom, depth_pin=depth_pin_nom, V_max=V, k_dielectric=k_dielectric_nom, t_dielectric=t_dielectric_nom,
#         #                                                               period=1, Fext=Fext)
#         T, Y, dY, X, dX, events, t_release_10pct, _, _, _, _, _, _ = sim_release_time(L_dielectric=L,
#                                                                                       w_pin=w_pin_nom, depth_pin=depth_pin_nom,
#                                                                                       V_max=V, t_dielectric=t_dielectric_nom,
#                                                                                       frequency=f_nom, Fext=Fext_nom,
#                                                                                       loadcell_start_ratio=loadcell_start_ratio,
#                                                                                       tau_1090=tau_1090)
#         t_engage.append(events)
#         # for Fext in Fext_range:
#         #     # T, Y, dY, t_end, Fk, Fb, Fes, Fconstant = sim_engagement_time(L_dielectric=L_dielectric_nom, w_pin=w_pin_nom, depth_pin=depth_pin_nom, V_max=V, k_dielectric=k_dielectric_nom, t_dielectric=t_dielectric_nom,
#         #     #                                                               period=1, Fext=Fext)
#         #     T, Y, dY, X, dX, events, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic = sim_engagement_time(L_dielectric=L_dielectric_nom,
#         #                                                                                                            w_pin=w_pin_nom, depth_pin=depth_pin_nom,
#         #                                                                                                            V_max=V, t_dielectric=t_dielectric_nom,
#         #                                                                                                            k_dielectric=k_dielectric_nom,
#         #                                                                                                            period=period_nom, Fext=Fext)
#         # t_engage.append(events)
#
#     # ax.plot(w_pin_range, t_engage, marker='o', label=r"$w_{pin}$" + " = {:0.1f} mm".format(w_pin * 1e3))
#     # ax4.semilogx(t_dielectric_range * 1e6, t_engage, label="{:d} V".format(int(V)))
#     ax4.semilogx(tau_1090_range * 1e6, t_engage, label="{:d} V".format(int(V)))
#     t_engage_Fext.append(t_engage)
#     # ax4.plot(Fext_range, t_engage, label="{:d} V".format(int(V)))
#     # ax4.plot(t_dielectric_range, t_engage, label="{:d} V".format(int(V)))
# # ax.set_xlabel("Voltage (V)")
# # ax4.set_xlabel(r"$T_d$ (μm)")
# # ax4.set_xlabel(r"$F_{preload}$ (N)")
# ax4.set_xlabel(r"Fall Time $t_{90\%-10\%}$ (µs)")
# ax4.set_xticks([1, 10, 1e2, 1e3, 1e4], [r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$'])
# fig.supylabel("Release Time (ms)", fontsize=20)
# # ax4.axvline(t_dielectric_nom * 1e6, c='k', ls='--')
# # ax.semilogy(True)
# ax4.grid(True)
# ax4.legend(fontsize=16)
# print("Results for tau1090:")
# print(tau_1090_range)
# print(t_engage_Fext)
# L_dielectric = 55.5e-3
# w_pin = 2e-3
# depth_pin = 2e-3
# V_max = 150
# k_dielectric = 54.2
# t_dielectric = 24e-6
# period = 1 / 1
# Fext = 0
# T, Y, dY, events, Fk, Fb, Fes, Fconstant = sim_engagement_time(L_dielectric=L_dielectric, w_pin=w_pin, depth_pin=depth_pin, V_max=V_max, k_dielectric=k_dielectric, t_dielectric=t_dielectric, period=period, Fext=Fext)
# print("Events:", events)
# fig, axs = plt.subplots(1, 2, layout='constrained', figsize=(11, 5))
# ax1, ax2 = axs
# ax1_right = ax1.twinx()
# ax1.plot(T * 1e6, Y * 1e6, '-o', c='tab:blue')
# ax1.set_xlabel("Time (us)")
# ax1.set_ylabel("Position (um)")
# # ax.axhline(t_air_final * 1e6, ls='--', c='k')
# ax1_right.plot(T * 1e6, dY, '-o', c='tab:orange')
# # ax_right.plot(sol.t * 1e6, [dy_dt(T, Y) for T, Y in zip(sol.t, np.transpose(sol.y))], '-o', c='tab:orange')
# ax1_right.set_ylabel("Velocity (m/s)")
# ax1.grid(True)
#
# ax2.plot(T * 1e6, Fes, label='Fes')
# ax2.plot(T * 1e6, Fk, label='Fk')
# ax2.plot(T * 1e6, -Fb, label='Fb')
# ax2.plot(T * 1e6, -Fconstant, label='Fconstant')
# ax2.set_xlabel("Time (us)")
# ax2.set_ylabel("Forces (N)")
# ax2.legend()
# ax2.grid(True)

print(timestamp)
np.save(save_folder + timestamp, [w_pin_range, t_engage_w_pin, t_dielectric_range, t_engage_t_dielectric, L_range, t_engage_L, tau_1090_range, t_engage_Fext], allow_pickle=True)

plt.savefig(save_folder + timestamp + ".png", dpi=300)
# # plt.savefig("figures/" + timestamp + ".svg")
plt.savefig(save_folder + timestamp + ".pdf")

end_time = datetime.now()
print("Total runtime:", (end_time - start_time).total_seconds(), 'seconds, or', (end_time - start_time).total_seconds() / 60, "minutes")
plt.show()
