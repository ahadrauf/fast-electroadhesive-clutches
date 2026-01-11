import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sim_engagement_time import *

now = datetime.now()
name_clarifier = "_engagement_time_sim_paramsweep"
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
fig, axs = plt.subplots(2, 2, layout='constrained', figsize=(6.4 * 1.2, 4.8 * 1.2))  # figsize=(10, 9))
ax1, ax2 = axs[0]
ax3, ax4 = axs[1]
L_dielectric_nom, w_pin_nom, depth_pin_nom, k_dielectric_nom, t_dielectric_nom, Fext_nom, f_nom = 55.5e-3, 2e-3, 2e-3, 54.2, 24e-6, 0.125, 1000  # F_extnom = 0.38
period_nom = 1 / 2 / f_nom
V_range = np.linspace(150, 300, 4)


def poly_fit(x, a, n, b):
    return a * np.power(x, n) + b


num_sim = 50
w_pin_range = np.power(10, np.linspace(-3, -1, num_sim))
t_engage_w_pin = []
for V in V_range:
    t_engage = []
    for w_pin in w_pin_range:
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

    ax1.semilogx(w_pin_range * 1e3, t_engage, label="{:d} V".format(int(V)))
    t_engage_w_pin.append(t_engage)

    popt, pcov = curve_fit(poly_fit, w_pin_range * 1e3, t_engage, p0=(1, 3, 0))
    print("Poly order:", popt[1], popt)
ax1.set_xlabel("Substrate Width (mm)")
ax1.axvline(2, c='k', ls='--')
ax1.axvline(2.5, c='k', ls='--')
ax1.axvline(3, c='k', ls='--')
ax1.axvline(4, c='k', ls='--')
ax1.axvline(5, c='k', ls='--')
ax1.axvline(6, c='k', ls='--')
ax1.set_xticks([1, 50, 100])
ax1.grid(True)

t_dielectric_range = np.power(10, np.linspace(-6, -2, num_sim))
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
        t_engage.append(events)
    ax2.semilogx(t_dielectric_range * 1e6, t_engage)
    t_engage_t_dielectric.append(t_engage)
ax2.set_xlabel("Dielectric Thickness (μm)")
ax2.axvline(t_dielectric_nom * 1e6, c='k', ls='--')
ax2.set_xticks([1, 10, 1e2, 1e3, 1e4], [r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$'])
ax2.tick_params(axis='x', which='minor', length=2, width=1, bottom=True)
ax2.grid(True)

L_range = np.power(10, np.linspace(-3, 0, num_sim))
t_engage_L = []
for V in V_range:
    t_engage = []
    for L in L_range:
        T, Y, dY, X, dX, events, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic = sim_engagement_time(L_dielectric=L,
                                                                                                               w_pin=w_pin_nom, depth_pin=depth_pin_nom,
                                                                                                               V_max=V, t_dielectric=t_dielectric_nom,
                                                                                                               k_dielectric=k_dielectric_nom,
                                                                                                               period=period_nom, Fext=Fext_nom,
                                                                                                               max_time_step=0.1e-6)
        t_engage.append(events)

    ax3.semilogx(L_range * 1e3, t_engage, label="{:d} V".format(int(V)))
    t_engage_L.append(t_engage)
ax3.set_xlabel(r"Substrate Length (mm)")
ax3.axvline(L_dielectric_nom * 1e3, c='k', ls='--')
ax3.set_xticks([1, 10, 1e2, 1e3], [r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$'])
ax3.grid(True)

tau_1090_range = np.power(10, np.linspace(-6, -2, num_sim))
t_engage_Fext = []
for V in V_range:
    t_engage = []
    for tau_1090 in tau_1090_range:
        T, Y, dY, X, dX, events, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic = sim_engagement_time(L_dielectric=L_dielectric_nom,
                                                                                                               w_pin=w_pin_nom, depth_pin=depth_pin_nom,
                                                                                                               V_max=V, t_dielectric=t_dielectric_nom,
                                                                                                               k_dielectric=k_dielectric_nom,
                                                                                                               period=period_nom, Fext=Fext_nom,
                                                                                                               max_time_step=0.1e-6, tau1090=tau_1090)
        t_engage.append(events)
    ax4.semilogx(tau_1090_range * 1e6, t_engage, label="{:d} V".format(int(V)))
    t_engage_Fext.append(t_engage)
ax4.set_xlabel(r"Rise Time $t_{10\%-90\%}$ (µs)")
ax4.set_xticks([1, 10, 1e2, 1e3, 1e4], [r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$'])
ax4.grid(True)
ax4.legend(fontsize=16)

fig.supylabel("Engagement Time (µs)", fontsize=20)

# np.save(save_folder + timestamp, [w_pin_range, t_engage_w_pin, t_dielectric_range, t_engage_t_dielectric, L_range, t_engage_L, tau_1090_range, t_engage_Fext])

end_time = datetime.now()
print("Total runtime:", (end_time - start_time).total_seconds(), 'seconds, or', (end_time - start_time).total_seconds() / 60, "minutes")
plt.show()
