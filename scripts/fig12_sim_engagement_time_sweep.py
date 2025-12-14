import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sim_engagement_time import *

now = datetime.now()
name_clarifier = "_engagement_time_sim_paramsweep"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
save_folder = "C:/Users/ahadrauf/Desktop/Research/latex/electroadhesive_dynamics_paper/figures_test/"
start_time = datetime.now()

plt.rcParams["font.family"] = "Arial"
plt.rc('font', size=24)  # controls default text size
plt.rc('axes', labelsize=20)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=20)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=20)  # fontsize of the y tick labels
plt.rc('legend', fontsize=20)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

file_loc = "../data/strain_tests/"

# Plot data
colors = ["#1964B0", "#F1932D", "#4DB264", "#DB060B"]  # "#F7F057" = yellow
fig, axs = plt.subplots(2, 2, layout='constrained', figsize=(6.4 * 1.2, 4.8 * 1.2))  # figsize=(10, 9))
ax1, ax2 = axs[0]
ax3, ax4 = axs[1]
L_dielectric_nom, w_pin_nom, depth_pin_nom, k_dielectric_nom, t_dielectric_nom, Fext_nom, f_nom = 55.5e-3, 2e-3, 2e-3, 54.2, 24e-6, 0.125, 1000  # F_extnom = 0.38
period_nom = 1 / 2 / f_nom
V_range = np.linspace(150, 300, 4)


# w_pin_range = np.array([1e-3, 2e-3, 3e-3, 4e-3, 5e-3, 6e-3, 7e-3])


def poly_fit(x, a, n, b):
    return a * np.power(x, n) + b


num_sim = 50
# w_pin_range = np.linspace(0.1e-3, 60e-3, num_sim)
w_pin_range = np.power(10, np.linspace(-3, -1, num_sim))
t_engage_w_pin = []
for V in V_range:
    t_engage = []
    for w_pin in w_pin_range:
        # T, Y, dY, t_end, Fk, Fb, Fes, Fconstant = sim_engagement_time(L_dielectric=L_dielectric_nom, w_pin=w_pin, depth_pin=depth_pin_nom, V_max=V, k_dielectric=k_dielectric_nom, t_dielectric=t_dielectric_nom,
        #                                                               period=1, Fext=Fext_nom)
        if w_pin < 1e-3:
            max_time_step = 0.1e-6
        elif w_pin > 30e-3:
            max_time_step = 10e-6
        else:
            max_time_step = 1e-6
        T, Y, dY, X, dX, events, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic = sim_engagement_time(L_dielectric=L_dielectric_nom,
                                                                                                               w_pin=w_pin, depth_pin=depth_pin_nom,
                                                                                                               V_max=V, t_dielectric=t_dielectric_nom,
                                                                                                               k_dielectric=k_dielectric_nom,
                                                                                                               period=period_nom, Fext=Fext_nom,
                                                                                                               max_time_step=max_time_step)
        t_engage.append(events)

    # ax.plot(w_pin_range, t_engage, marker='o', label=r"$w_{pin}$" + " = {:0.1f} mm".format(w_pin * 1e3))
    # ax1.plot(w_pin_range * 1e3, t_engage, label="{:d} V".format(int(V)))
    ax1.semilogx(w_pin_range * 1e3, t_engage, label="{:d} V".format(int(V)))
    t_engage_w_pin.append(t_engage)

    popt, pcov = curve_fit(poly_fit, w_pin_range * 1e3, t_engage, p0=(1, 3, 0))
    print("Poly order:", popt[1], popt)
# ax.set_xlabel("Voltage (V)")
ax1.set_xlabel(r"$w_s$ (mm)")
# ax1.set_ylabel("Engagement Time (µs)")
ax1.axvline(2, c='k', ls='--')
ax1.axvline(2.5, c='k', ls='--')
ax1.axvline(3, c='k', ls='--')
ax1.axvline(4, c='k', ls='--')
ax1.axvline(5, c='k', ls='--')
ax1.axvline(6, c='k', ls='--')
# ax1.set_xticks([0, 20, 40, 60])
ax1.set_xticks([1, 10, 1e2], [r'$10^0$', r'$10^1$', r'$10^2$'])
# ax.semilogy(True)
ax1.grid(True)
# ax1.legend()

# k_range = np.linspace(3, 100, num_sim)
# f_range = np.power(10, np.linspace(0, 5, 2 * num_sim))
t_dielectric_range = np.power(10, np.linspace(-6, -2, num_sim))  # 10000e-6
t_engage_t_dielectric = []
for V in V_range:
    t_engage = []
    for t_dielectric in t_dielectric_range:
        T, Y, dY, X, dX, events, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic = sim_engagement_time(L_dielectric=L_dielectric_nom,
                                                                                                               w_pin=w_pin_nom, depth_pin=depth_pin_nom,
                                                                                                               V_max=V, t_dielectric=t_dielectric,
                                                                                                               k_dielectric=k_dielectric_nom,
                                                                                                               period=period_nom, Fext=Fext_nom,
                                                                                                               max_time_step=0.25e-6)
        # for k in k_range:
        #     T, Y, dY, X, dX, events, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic = sim_engagement_time(L_dielectric=L_dielectric_nom,
        #                                                                                                            w_pin=w_pin_nom, depth_pin=depth_pin_nom,
        #                                                                                                            V_max=V, t_dielectric=t_dielectric_nom,
        #                                                                                                            k_dielectric=k,
        #                                                                                                            period=1 / 2 / f_nom, Fext=Fext_nom)
        #     t_engage.append(events)
        # for f in f_range:
        #     # T, Y, dY, t_end, Fk, Fb, Fes, Fconstant = sim_engagement_time(L_dielectric=L_dielectric_nom, w_pin=w_pin_nom, depth_pin=depth_pin_nom, V_max=V, k_dielectric=k_dielectric_nom, t_dielectric=t_dielectric_nom,
        #     #                                                               period=1 / f, Fext=Fext_nom)
        #     T, Y, dY, X, dX, events, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic = sim_engagement_time(L_dielectric=L_dielectric_nom,
        #                                                                                                            w_pin=w_pin_nom, depth_pin=depth_pin_nom,
        #                                                                                                            V_max=V, t_dielectric=t_dielectric_nom,
        #                                                                                                            k_dielectric=k_dielectric_nom,
        #                                                                                                            period=1 / 2 / f, Fext=Fext_nom)
        t_engage.append(events)

    # ax.plot(w_pin_range, t_engage, marker='o', label=r"$w_{pin}$" + " = {:0.1f} mm".format(w_pin * 1e3))
    # ax2.plot(k_range, t_engage, label="{:d} V".format(int(V)))
    ax2.semilogx(t_dielectric_range * 1e6, t_engage)
    # ax2.semilogx(f_range, t_engage, label="{:d} V".format(int(V)))

    t_engage_t_dielectric.append(t_engage)
# ax.set_xlabel("Voltage (V)")
# ax2.set_xlabel(r"$\kappa_s$ Dielectric Constant")
ax2.set_xlabel(r"$T_d$ (μm)")
ax2.axvline(t_dielectric_nom * 1e6, c='k', ls='--')
# ax2.axvline(k_dielectric_nom, c='k', ls='--')
# ax2.set_xlabel("Drive Frequency (Hz)")
# ax2.set_ylabel("Engagement Time (µs)")
# for f in [1, 10, 100, 1000, 2000, 5000, 10000]:
#     ax2.axvline(f, c='k', ls='--')
# # ax.semilogy(True)
# ax2.set_xticks([1, 10, 1e2, 1e3, 1e4, 1e5], [r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$'])
ax2.set_xticks([1, 10, 1e2, 1e3, 1e4], [r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$'])
# ax2.xaxis.set_major_formatter(mtick.ScalarFormatter(useMathText=True))
# ax2.xaxis.set_minor_locator(plt.LogLocator(subs=np.arange(2, 10)))
ax2.tick_params(axis='x', which='minor', length=2, width=1, bottom=True)
ax2.grid(True)
# ax2.legend()

L_range = np.power(10, np.linspace(-3, 0, num_sim))
t_engage_L = []
for V in V_range:
    t_engage = []
    for L in L_range:
        # T, Y, dY, t_end, Fk, Fb, Fes, Fconstant = sim_engagement_time(L_dielectric=L, w_pin=w_pin_nom, depth_pin=depth_pin_nom, V_max=V, k_dielectric=k_dielectric_nom, t_dielectric=t_dielectric_nom,
        #                                                               period=1, Fext=Fext_nom)
        T, Y, dY, X, dX, events, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic = sim_engagement_time(L_dielectric=L,
                                                                                                               w_pin=w_pin_nom, depth_pin=depth_pin_nom,
                                                                                                               V_max=V, t_dielectric=t_dielectric_nom,
                                                                                                               k_dielectric=k_dielectric_nom,
                                                                                                               period=period_nom, Fext=Fext_nom,
                                                                                                               max_time_step=0.1e-6)
        t_engage.append(events)

    # ax.plot(w_pin_range, t_engage, marker='o', label=r"$w_{pin}$" + " = {:0.1f} mm".format(w_pin * 1e3))
    ax3.semilogx(L_range * 1e3, t_engage, label="{:d} V".format(int(V)))
    t_engage_L.append(t_engage)
# ax.set_xlabel("Voltage (V)")
ax3.set_xlabel(r"$L_s$ (mm)")
# ax3.set_ylabel("Engagement Time (µs)")
ax3.axvline(L_dielectric_nom * 1e3, c='k', ls='--')
ax3.set_xticks([1, 10, 1e2, 1e3], [r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$'])
# ax.semilogy(True)
ax3.grid(True)
# ax3.legend()

# t_dielectric_range = np.linspace(1e-6, 1000e-6, num_sim)  # 10000e-6
# t_dielectric_range = np.power(10, np.linspace(-6, -2, num_sim))  # 10000e-6
# Fext_range = np.linspace(0, 1, num_sim)
tau_1090_range = np.power(10, np.linspace(-6, -2, num_sim))
t_engage_Fext = []
for V in V_range:
    t_engage = []
    # for t_dielectric in t_dielectric_range:
    #     T, Y, dY, X, dX, events, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic = sim_engagement_time(L_dielectric=L_dielectric_nom,
    #                                                                                                            w_pin=w_pin_nom, depth_pin=depth_pin_nom,
    #                                                                                                            V_max=V, t_dielectric=t_dielectric,
    #                                                                                                            k_dielectric=k_dielectric_nom,
    #                                                                                                            period=period_nom, Fext=Fext_nom)
    # for Fext in Fext_range:
    #     # T, Y, dY, t_end, Fk, Fb, Fes, Fconstant = sim_engagement_time(L_dielectric=L_dielectric_nom, w_pin=w_pin_nom, depth_pin=depth_pin_nom, V_max=V, k_dielectric=k_dielectric_nom, t_dielectric=t_dielectric_nom,
    #     #                                                               period=1, Fext=Fext)
    #     T, Y, dY, X, dX, events, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic = sim_engagement_time(L_dielectric=L_dielectric_nom,
    #                                                                                                            w_pin=w_pin_nom, depth_pin=depth_pin_nom,
    #                                                                                                            V_max=V, t_dielectric=t_dielectric_nom,
    #                                                                                                            k_dielectric=k_dielectric_nom,
    #                                                                                                            period=period_nom, Fext=Fext,
    #                                                                                                            max_time_step=0.1e-6)
    for tau_1090 in tau_1090_range:
        # T, Y, dY, t_end, Fk, Fb, Fes, Fconstant = sim_engagement_time(L_dielectric=L_dielectric_nom, w_pin=w_pin_nom, depth_pin=depth_pin_nom, V_max=V, k_dielectric=k_dielectric_nom, t_dielectric=t_dielectric_nom,
        #                                                               period=1, Fext=Fext)
        T, Y, dY, X, dX, events, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic = sim_engagement_time(L_dielectric=L_dielectric_nom,
                                                                                                               w_pin=w_pin_nom, depth_pin=depth_pin_nom,
                                                                                                               V_max=V, t_dielectric=t_dielectric_nom,
                                                                                                               k_dielectric=k_dielectric_nom,
                                                                                                               period=period_nom, Fext=Fext_nom,
                                                                                                               max_time_step=0.1e-6, tau1090=tau_1090)
        t_engage.append(events)

    # ax.plot(w_pin_range, t_engage, marker='o', label=r"$w_{pin}$" + " = {:0.1f} mm".format(w_pin * 1e3))
    # ax4.semilogx(t_dielectric_range * 1e6, t_engage, label="{:d} V".format(int(V)))
    # ax4.plot(Fext_range, t_engage, label="{:d} V".format(int(V)))
    ax4.semilogx(tau_1090_range * 1e6, t_engage, label="{:d} V".format(int(V)))
    t_engage_Fext.append(t_engage)
# ax.set_xlabel("Voltage (V)")
# ax4.set_xlabel(r"$t_{substrate}$ (μm)")
# ax4.set_xlabel(r"$F_{preload}$ (N)")
ax4.set_xlabel(r"Rise Time $t_{10\%-90\%}$ (µs)")
ax4.set_xticks([1, 10, 1e2, 1e3, 1e4, 1e5], [r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$'])
fig.supylabel("Engagement Time (µs)", fontsize=20)
# ax4.axvline(t_dielectric_nom * 1e6, c='k', ls='--')
# ax.semilogy(True)
ax4.grid(True)
ax4.legend(fontsize=16)

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

# plt.savefig(save_folder + timestamp + ".png", dpi=300)
# # # plt.savefig("figures/" + timestamp + ".svg")
# plt.savefig(save_folder + timestamp + ".pdf")
#
# # np.save(save_folder + timestamp, [w_pin_range, t_engage_w_pin, t_dielectric_range, t_engage_t_dielectric, L_range, t_engage_L, Fext_range, t_engage_Fext])
# np.save(save_folder + timestamp, [w_pin_range, t_engage_w_pin, t_dielectric_range, t_engage_t_dielectric, L_range, t_engage_L, tau_1090_range, t_engage_Fext])

end_time = datetime.now()
print("Total runtime:", (end_time - start_time).total_seconds(), 'seconds, or', (end_time - start_time).total_seconds() / 60, "minutes")
plt.show()
