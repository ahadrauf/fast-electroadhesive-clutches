import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from scipy.integrate import solve_ivp, quad
from scipy.special import kv
from scipy.optimize import root_scalar
from datetime import datetime
import csv

now = datetime.now()
name_clarifier = "_engagement_time_sim_paramsweep"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
save_folder = "C:/Users/ahadrauf/Desktop/Research/latex/electroadhesive_dynamics_paper/figures_test/"
start_time = datetime.now()

plt.rcParams["font.family"] = "Arial"
# plt.rcParams["font.family"] = "serif"
# plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rc('font', size=16)  # controls default text size
plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=16)  # fontsize of the y tick labels
plt.rc('legend', fontsize=16)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
# plt.style.use('tableau-colorblind10')

N = 28
eps0 = 8.854e-12
# k_dielectric = 54.2  # 50
w_trace = 1.5e-3
s_trace = 0.5e-3
# t_dielectric = 24e-6
# w_dielectric = 2e-3  # 2e-3
L_pin = 100e-3
g = 9.81  # m/s^2
E = 150e6
# V_max = 150
# print("Mass pin:", m_pin, "Mass PVDF:", m_pvdf)

mu_objet_s = 0.281  # 0.173  # 0.28
mu_pvdf_s = 0.188  # 0.154  # 0.38
mu_objet_k = 0.173
mu_pvdf_k = 0.154
mu_tot_k = mu_objet_k + mu_pvdf_k
mu_tot_s = mu_objet_s + mu_pvdf_s
t_air_initial = 20e-6  # 38e-6  # 45e-6  # 34.2e-6  # t_air_contact - (m_pvdf * g + F_preload) / k_contact
k_lc = 18.017806799339144 * 1e3  # N/mm
omega_lc = 2 * np.pi * 643.4036  # calcualted on 20241029 Daily Notes
m_lc = k_lc / np.square(omega_lc)  # N/mm --> N/m / (1/s^2)
# print("Load cell effective mass", m_lc, "kg, Expected quarter period", 2 * np.pi / omega_lc / 4 * 1e6, "us")

with open('C:/Users/ahadrauf/OneDrive - Stanford/Research/Papers/TMech 2024/cole_cole_inverse_laplace_pvdftrfecfe.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    t_cc = []
    e_cc = []
    for row in reader:
        t_cc.append(float(row[0]))
        e_cc.append(float(row[4]))
    t_cc = np.array(t_cc)
    e_cc = np.array(e_cc)
e_cc_max = np.max(e_cc)

h2over4 = lambda h: np.square(h) / 4
h2 = lambda h: np.square(h)
F32 = lambda h: np.exp(-h2over4(h)) * np.sqrt(h) * ((h2(h) + 1) * kv(0.25, h2over4(h)) - h2(h) * kv(0.75, h2over4(h))) / 4 / np.sqrt(np.pi)
# sigma_gw, C_gw = 2.38436937e-06, 2224250654027140.2  # 2.40459756e+15  # 3.13725051e+11
# sigma_gw, C_gw = 2.80263157e-06, 5.39050698e+14
sigma_gw, C_gw = 2.01077514e-06, 3.26077705e+15  # 2mm, all data range
# sigma_gw, C_gw = 2.369088099456535e-06, 2086492814152674.0  # avg, all data range
# sigma_gw, C_gw = 2.76314774e-06, 3.05500533e+14
Fes_scalar = 4.0163678631132145


def sim_engagement_time(L_dielectric, w_pin, depth_pin, V_max, k_dielectric, t_dielectric, period, Fext, max_time_step=1e-6, tau1090=8.2e-6):
    w_dielectric = w_pin
    m_pin = 8.049e3 * w_pin * depth_pin * L_pin
    m_pvdf = 1.78e3 * t_dielectric * L_dielectric * w_dielectric  # 1.78 g/cc
    m_tot = m_pin + m_pvdf
    # t_air_contact = 10e-6
    Neff = N * (L_dielectric / 55.5e-3)
    # t_air_initial_after_weight = t_air_initial  # - np.sqrt(2 * (m_pvdf * g + Fext) * (t_air_initial - t_air_contact) / k_contact)
    # print("Actual t_air_initial", 1e6 * t_air_initial_after_weight, "vs. original", t_air_initial * 1e6)
    # print("Contact at", t_air_contact * 1e6, "um")
    F_gw = lambda h: C_gw * np.power(sigma_gw, 1.5) * w_pin * L_dielectric * F32(h / sigma_gw)
    # print("F_gw at 5 um", F_gw(5e-6))
    # Fes_scalar = 1.70824E-05 * (V_max**2) - 0.01352372 * V_max + 3.208127237
    F_es = lambda k, t_air, V: Fes_scalar / 8 * eps0 * np.square(k) * (V**2) * (Neff * w_trace * w_pin) / np.square(t_dielectric + k * t_air)
    # F_es = lambda k, t_air, V: Fes_scalar / 8 * eps0 * np.square(k) * (V**2) * (N * w_trace * w_pin) / np.sqrt(2 * np.pi) / sigma_gw * \
    #                            quad(lambda y: np.divide(np.exp(-np.square(t_air - y) / 2 / np.square(sigma_gw)), np.square(t_dielectric + k * y)), -t_dielectric / k * 0.83, 1e-3)[0]
    V_t = lambda t: V_max * (1 - np.exp(-t / (tau1090 / np.log(9))))  # conversion from 10-90 rise time to tau
    epsr = lambda t: np.interp(t, t_cc, e_cc)

    t_air_initial_after_weight = root_scalar(lambda t_air: F_gw(t_air) - m_pvdf * g - Fext, x0=5e-6, x1=4e-6).root
    t_air_final = root_scalar(lambda t_air: F_gw(t_air) - F_es(k=e_cc_max, t_air=t_air, V=V_max) - m_pvdf * g - Fext,
                              x0=5e-6, x1=4e-6).root
    print("Initial air gap:", t_air_initial_after_weight * 1e6, "um")
    print("Predicted final air gap:", t_air_final * 1e6, "um")

    made_contact = lambda t, Y: Y[0] - (t_air_initial_after_weight + 0.005 * (t_air_final - t_air_initial_after_weight))  # 0.9
    made_contact.terminal = True
    made_contact.direction = 0

    def calc_forces(t_curr, y, ydot):
        # Fk = 0 * y
        Fk = F_gw(y)
        # Fk = k_contact / (t_air_initial - t_air_contact) * 0.5 * np.square(t_air_initial - y)
        S1, S2 = min(w_dielectric, L_dielectric), max(w_dielectric, L_dielectric)
        b_sf = 96 / (np.pi**4) * 1.82e-5 * S2 * (S1**3) / np.power(y, 3) * (1 - 0.58 * (S1 / S2))
        Fb = b_sf * ydot
        k_t = epsr(t_curr % period)

        if type(y) == np.float64:
            Fes = F_es(k=k_t, t_air=y, V=V_t(t_curr))
        else:
            Fes = [F_es(k=k_ti, t_air=yi, V=V_t(t_curri)) for k_ti, yi, t_curri in zip(k_t, y, t_curr)]
        Fconstant = m_pvdf * g * np.ones_like(y) + Fext
        Nb = m_pin * g + m_pvdf * g + Fext
        Nf = Fk
        Flinear_static = mu_objet_s * Nb + mu_pvdf_s * Nf
        Flinear_kinetic = mu_objet_k * Nb * np.ones_like(y)  # + mu_pvdf_k * Nf
        return Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic

    def dy_dt(t, Y):
        y, ydot, x, xdot = Y  # x and xdot are included for compatibility with sim_release_time (they're not actually used in engagement_time sims)
        Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic = calc_forces(t, y, ydot)
        yddot = 1 / m_pvdf * (Fk - Fb - Fes - Fconstant)
        xddot = 0
        return [ydot, yddot, xdot, xddot]

    # if t_air_initial_after_weight > 100e-6:
    #     print("Original initial air gap after weight", t_air_initial_after_weight)
    #     t_air_initial_after_weight = 10e-6

    Y0 = [t_air_initial_after_weight, 0, 0, 0]
    sol = solve_ivp(dy_dt, [0.1e-8, 1000e-6], Y0, max_step=max_time_step, events=[made_contact])

    T, Y, dY, X, dX = sol.t, sol.y[0, :], sol.y[1, :], sol.y[2, :], sol.y[3, :]
    Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic = calc_forces(T, Y, dY)
    print("Final gap: {:0.3e}, Fes = {:0.3e}, Fk = {:0.3e}, Fb = {:0.3e}, Fconstant = {:0.3e}".format(Y[-1], Fes[-1], Fk[-1], Fb[-1], Fconstant[-1]))

    if len(sol.t_events[0]) > 0:
        print("L = {:0.1f}mm, w_pin = {:0.1f}mm, V = {:0.1f}, k = {:0.1f}, t_d = {:0.1f}um --> t_engage = {:0.3f} us".format(L_dielectric * 1e3, w_pin * 1e3, V_max, k_dielectric, t_dielectric * 1e6, sol.t_events[0][0] * 1e6))
        return T, Y, dY, X, dX, sol.t_events[0][0] * 1e6, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic
    else:
        print("L = {:0.1f}mm, w_pin = {:0.1f}mm, V = {:0.1f}, k = {:0.1f}, t_d = {:0.1f}um --> t_engage = Ran until endstop".format(L_dielectric * 1e3, w_pin * 1e3, V_max, k_dielectric, t_dielectric * 1e6))
        return T, Y, dY, X, dX, np.nan, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic
    # print(sol.y)


if __name__ == '__main__':
    L_dielectric = 55.5e-3
    w_pin = 2e-3
    depth_pin = 2e-3
    V_max = 300
    k_dielectric = 54.2
    t_dielectric = 24e-6
    frequency = 1000
    period = 1 / 2 / frequency
    Fext = 0  # 0.04273761595304056
    T, Y, dY, X, dX, events, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic = sim_engagement_time(L_dielectric=L_dielectric, w_pin=w_pin,
                                                                                                           depth_pin=depth_pin, V_max=V_max,
                                                                                                           t_dielectric=t_dielectric, period=period,
                                                                                                           Fext=Fext, k_dielectric=k_dielectric)
    print("Events:", events)
    fig, axs = plt.subplots(1, 2, layout='constrained', figsize=(11, 5))
    ax1, ax2 = axs
    ax1_right = ax1.twinx()
    ax1.plot(T * 1e6, Y * 1e6, '-o', c='tab:blue', label=r"$y$")
    ax1.plot(T * 1e6, X * 1e6, '-o', c='tab:green', label=r"$x$")
    ax1.axhline(X[0] * 1e6, c='k', ls='--')
    ax1.set_xlabel("Time (us)")
    ax1.set_ylabel("Position (um)")
    # ax.axhline(t_air_final * 1e6, ls='--', c='k')
    ax1_right.plot(T * 1e6, dY * 1e3, '-o', c='tab:orange', label=r"$\dot{y}$")
    ax1_right.plot(T * 1e6, dX * 1e3, '-o', c='tab:brown', label=r"$\dot{x}$")
    ax1.legend(loc='upper right')
    ax1_right.legend(loc='lower right')
    ax1_right.set_ylabel("Velocity (mm/s)")
    ax1.grid(True)

    ax2.plot(T * 1e6, Fes, label='Fes')
    ax2.plot(T * 1e6, Fk, label='Fk')
    ax2.plot(T * 1e6, Fb, label='Fb')
    ax2.plot(T * 1e6, Fconstant, label='Fconstant')
    ax2.plot(T * 1e6, k_lc * X, label='Load Cell Spring')
    ax2.plot(T * 1e6, Flinear_static, label='Ffriction_s')
    ax2.plot(T * 1e6, Flinear_kinetic, label='Ffriction_k')
    ax2.set_xlabel("Time (us)")
    ax2.set_ylabel("Forces (N)")
    ax2.legend()
    ax2.grid(True)

    plt.show()
