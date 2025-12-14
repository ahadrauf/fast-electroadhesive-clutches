import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, root_scalar
from scipy.integrate import solve_ivp, quad
from scipy.special import kv
import csv

N = 28
eps0 = 8.854e-12
w_trace = 1.5e-3
s_trace = 0.5e-3
L_pin = 95e-3
g = 9.81  # m/s^2

mu_objet_s = 0.281  # 0.173  # 0.28
mu_pvdf_s = 0.188  # 0.154  # 0.38
mu_objet_k = 0.173
mu_pvdf_k = 0.154
mu_tot_k = mu_objet_k + mu_pvdf_k
mu_tot_s = mu_objet_s + mu_pvdf_s
t_air_initial = 10e-6  # 38e-6  # 45e-6  # 34.2e-6  # t_air_contact - (m_pvdf * g + F_preload) / k_contact
k_lc = 18.017806799339144 * 1e3  # N/mm
omega_lc = 2 * np.pi * 643.4036  # calcualted on 20241029 Daily Notes
m_lc = 0.001100401067840215  # k_lc / np.square(omega_lc)  # N/mm --> N/m / (1/s^2)
b_lc = 0.42743335268958244  # 0.3860206506737978

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
# sigma_gw, C_gw = 2.38436937e-06, 2224250654027140.2  # 1183770965076373.5  # 2224250654027140.2  # 2.40459756e+15  # 3.13725051e+11  # 2.40459756e+12  #
sigma_gw, C_gw = 2.80263157e-06, 5.39050698e+14


def sim_release_time(L_dielectric, w_pin, depth_pin, V_max, t_dielectric, frequency, Fext, loadcell_start_ratio=1, tau_1090=5.2e-6, Fes_scalar_ideal=False):
    w_dielectric = w_pin
    m_pin = 8.049e3 * w_pin * depth_pin * L_pin
    m_pvdf = 1.78e3 * t_dielectric * L_dielectric * w_dielectric  # 1.78 g/cc
    Neff = N * (L_dielectric / 55.5e-3)
    if Fes_scalar_ideal:
        Fes_scalar = 1
    else:
        Fes_scalar = 1.70824E-05 * (V_max**2) - 0.01352372 * V_max + 3.208127237
    F_es = lambda k, t_air, V: Fes_scalar / 8 * eps0 * np.square(k) * (V**2) * (Neff * w_trace * w_pin) / np.sqrt(2 * np.pi) / sigma_gw * \
                               quad(lambda y: np.divide(np.exp(-np.square(t_air - y) / 2 / np.square(sigma_gw)), np.square(t_dielectric + k * y)), -t_dielectric / k * 0.83, 1e-3)[0]
    F_gw = lambda h: C_gw * np.power(sigma_gw, 1.5) * (Neff * w_trace * w_pin) * F32(h / sigma_gw)  # L_dielectric

    V_t = lambda t: V_max * np.exp(-t / (tau_1090 / np.log(9)))  # conversion from 10-90 rise time to tau
    epsr = lambda t: e_cc_max - np.interp(t, t_cc, e_cc)
    e_high, e_low, tau, alpha = 4, 54.2, 2.82e-6, 0.438
    cole_cole = lambda f: np.real(e_high + (e_low - e_high) / (1 + np.power(1j * 2 * np.pi * f * tau, 1 - alpha)))
    eps_initial = cole_cole(frequency)

    t_air_initial_after_weight = root_scalar(lambda t_air: F_gw(t_air) - F_es(k=eps_initial, t_air=t_air, V=V_max) - m_pvdf * g - Fext, x0=4e-6, x1=5e-6).root
    t_air_final_expected = root_scalar(lambda t_air: F_gw(t_air) - m_pvdf * g - Fext, x0=4e-6, x1=5e-6).root

    Fes_initial = F_es(k=eps_initial, t_air=t_air_initial_after_weight, V=V_max)
    # print("Load cell initial start ratio", loadcell_start_ratio)
    Flinear_static_initial = mu_objet_s * (m_pin * g + m_pvdf * g + Fext) + mu_pvdf_s * F_gw(t_air_initial_after_weight)
    Flinear_static_final_expected = mu_objet_s * (m_pin * g + m_pvdf * g + Fext) + mu_pvdf_s * F_gw(t_air_final_expected)
    loadcell_dx_t0 = Flinear_static_initial / k_lc  # mm
    loadcell_dx_final = Flinear_static_final_expected / k_lc
    loadcell_dx_initial = loadcell_dx_final + (loadcell_dx_t0 - loadcell_dx_final) * (1 - loadcell_start_ratio)  # Flinear_static_initial / Flinear_static_initial_ideal
    Y0 = [t_air_initial_after_weight, 0, loadcell_dx_initial, 0]

    # print("Actual t_air_initial", 1e6*t_air_initial_after_weight, "final", t_air_final_expected*1e6)
    # print("Initial load cell dx", loadcell_dx_t0*1e6, "final", loadcell_dx_final*1e6, "started at", loadcell_dx_initial*1e6)
    # print("Predicted max force:", Flinear_static_initial - Flinear_static_final_expected, "N", (Flinear_static_initial - Flinear_static_final_expected)/w_pin/60, "kPa")

    def made_contact(t, Y):
        if (Y[3] > -1e-5) and (Y[2] < loadcell_dx_initial - 1e-8):  # (t > 100e-6):
            return -1
        else:
            return 1

    made_contact.terminal = True
    made_contact.direction = -1

    def calc_forces(t_curr, y, ydot):
        Fk = F_gw(y)
        S1, S2 = min(w_dielectric, L_dielectric), max(w_dielectric, L_dielectric)
        b_sf = 96 / (np.pi**4) * 1.82e-5 * S2 * (S1**3) / np.power(y, 3) * (1 - 0.58 * (S1 / S2))
        Fb = -b_sf * ydot  # < 0 (because ydot > 0)
        if type(y) == np.float64:
            Fes = F_es(k=epsr(t_curr), t_air=y, V=V_t(t_curr))
        else:
            Fes = [F_es(k=epsr(t_curr_i), t_air=y_i, V=V_t(t_curr_i)) for t_curr_i, y_i in zip(t_curr, y)]
        Fconstant = m_pvdf * g * np.ones_like(y) + Fext
        yddot = 1 / m_pvdf * (Fk + Fb - Fes - Fconstant)
        Nb = m_pin * g + m_pvdf * g + Fext + m_pvdf * yddot  # - Fb
        Nf = Fk
        Flinear_static = mu_objet_s * Nb + mu_pvdf_s * Nf
        Flinear_kinetic = mu_objet_k * Nb + mu_pvdf_k * Nf
        return Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic

    def dy_dt(t, Y):
        y, ydot, x, xdot = Y
        Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic = calc_forces(t, y, ydot)
        yddot = 1 / m_pvdf * (Fk + Fb - Fes - Fconstant)
        if abs(xdot) < 1e-8 and k_lc * x + b_lc * xdot <= Flinear_static:
            xddot = 0
        else:
            xddot = 1 / m_lc * (-Flinear_kinetic * np.sign(xdot) - k_lc * x - b_lc * xdot)
        return [ydot, yddot, xdot, xddot]

    # print("Initial Fes", Fes_initial, "Linear static friction", Flinear_static_initial, "Loadcell dx", loadcell_dx_initial*1e6, "um")
    # print("Initial position:", Y0, Fext)
    sol = solve_ivp(dy_dt, [0.1e-6, 10e-3], Y0, max_step=1e-6, events=[made_contact])

    T, Y, dY, X, dX = sol.t, sol.y[0, :], sol.y[1, :], sol.y[2, :], sol.y[3, :]
    Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic = calc_forces(T, Y, dY)

    threshold_10pct = X[-1] + (X[0] - X[-1]) * 0.9
    release_time_10pct = T[np.argmin(np.abs(X - threshold_10pct))]
    # print("Final gap: {:0.3e}, Final Load Cell x: {:0.3e}".format(Y[-1], X[-1]))

    # Can delete this for most of code
    # t_air_with_Fes = root_scalar(lambda t_air: F_gw(t_air) - F_es(k=e_cc_max, t_air=t_air, V=V_max) - m_pvdf * g - Fext, bracket=(1e-6, 10e-6)).root
    # t_air_without_Fes = root_scalar(lambda t_air: F_gw(t_air) - m_pvdf * g - Fext, x0=3e-6, x1=4e-6).root
    # Fshear_pred = (mu_pvdf_s * F_gw(t_air_with_Fes) - mu_pvdf_k * F_gw(t_air_without_Fes))  #  / (w_pin * L_dielectric) / 1e3

    if len(sol.t_events[0]) > 0:
        print("L = {:0.3f}mm, w_pin = {:0.3f}mm, V = {:0.1f}, frequency = {:0.1f}, t_d = {:0.1f}um, F_preload = {:0.3f}um --> t_release = {:0.6f} us".format(L_dielectric * 1e3, w_pin * 1e3, V_max,
                                                                                                                                                             frequency, t_dielectric * 1e6, Fext,
                                                                                                                                                             sol.t_events[0][0] * 1e6))
        # print("L = {:0.3f}mm, w_pin = {:0.3f}mm, V = {:0.1f}, frequency = {:0.1f}, t_d = {:0.1f}um, F_preload = {:0.3f} --> t_release = {:0.6f} us, Fshear = {:0.6f} N".format(L_dielectric * 1e3, w_pin * 1e3, V_max,
        #                                                                                                                                                                        frequency, t_dielectric * 1e6, Fext,
        #                                                                                                                                                                        sol.t_events[0][0] * 1e6, Fshear_pred))
        # return T, Y, dY, X, dX, sol.t_events[0][0] * 1e3, tr10pct, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic
        return T, Y, dY, X, dX, sol.t_events[0][0] * 1e3, release_time_10pct * 1e3, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic
        # return T, Y, dY, X, dX, sol.t_events[0][0]*1e3, release_time_10pct*1e3, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    else:
        print("L = {:0.1f}mm, w_pin = {:0.1f}mm, V = {:0.1f}, frequency = {:0.1f}, t_d = {:0.1f}um, F_preload = {:0.1f}um --> t_release = Ran until endstop".format(L_dielectric * 1e3, w_pin * 1e3, V_max,
                                                                                                                                                                    frequency, t_dielectric * 1e6, Fext))
        # return T, Y, dY, X, dX, np.nan, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic
        # return T, Y, dY, X, dX, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        return T, Y, dY, X, dX, 1000, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    # print(sol.y)


if __name__ == '__main__':
    L_dielectric = 55.5e-3
    w_pin = 2e-3
    depth_pin = 2e-3
    V_max = 300
    k_dielectric = 54.2
    t_dielectric = 24e-6
    frequency = 1000
    load_cell_start_ratio = 0.8
    Fext = 1  # 0.08439257560224675  # 0.14077977424369925  # 0.09368685112003708  # 0.04273761595304056
    T, Y, dY, X, dX, events, tr10pct, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic = sim_release_time(L_dielectric=L_dielectric, w_pin=w_pin, depth_pin=depth_pin, V_max=V_max, t_dielectric=t_dielectric, frequency=frequency, Fext=Fext, loadcell_start_ratio=load_cell_start_ratio)
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
    ax1_right.plot(T * 1e6, dY, '-o', c='tab:orange', label=r"$\dot{y}$")
    ax1_right.plot(T * 1e6, dX, '-o', c='tab:brown', label=r"$\dot{x}$")
    ax1.legend(loc='upper right')
    ax1_right.legend(loc='lower right')
    # ax_right.plot(sol.t * 1e6, [dy_dt(T, Y) for T, Y in zip(sol.t, np.transpose(sol.y))], '-o', c='tab:orange')
    ax1_right.set_ylabel("Velocity (m/s)")
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
