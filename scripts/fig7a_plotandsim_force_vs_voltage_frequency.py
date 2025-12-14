import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit, root_scalar
from scipy.integrate import quad
from scipy.special import kv
from datetime import datetime

now = datetime.now()
name_clarifier = "_force_capacity_vs_voltage_frequency"
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

file_name = "20241114_14_14_33_force_capacity_vs_voltage_frequency_newwheatstonedata_plotdata"

data = np.load("../data/" + file_name + ".npy", allow_pickle=True)
all_V, all_freq, all_x1, all_y1, all_yerr1, all_x2, all_y2, all_yerr2, all_load_cell_disengage_forces, all_load_cell_preload_forces, all_load_cell_max_force = data

fig2, ax2 = plt.subplots(1, 1, figsize=(6, 4), layout='constrained')
labels = ["DC", r"$10^0$ Hz", r"$10^1$ Hz", r"$10^2$ Hz", r"$10^3$ Hz", r"$10^4$ Hz", r"$10^4$ Hz (Model)"]
ax1_lines = []
markers = ['o', '^', 'd', 'x']
colors = ["#1964B0", "#F1932D", "#4DB264", "#DB060B"]


def exponential_func(v, a, n):
    return a*np.power(v, n)


def fit_to_curve(v, f):
    if v[0] != 0:
        v = np.insert(v, 0, 0)
        f = np.insert(f, 0, 0)
    popt, pcov = curve_fit(exponential_func, v, f, bounds=([0, 0.5], [100, 4]))
    ss_res = np.sum(np.square(f - exponential_func(v, *popt)))
    ss_tot = np.sum(np.square(f - np.mean(f)))
    r2 = r2_score(exponential_func(v, *popt), f)  # np.corrcoef(f, exponential_func(v, *popt))[0, 1]**2
    r2 = 1 - ss_res/ss_tot
    # mape = 100/len(v) * np.sum(np.abs(np.divide(f - exponential_func(v, *popt), f)))
    return popt, r2


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

    mu_objet_s = 0.281  # 0.173  # 0.28
    mu_pvdf_s = 0.188  # 0.154  # 0.38
    mu_objet_k = 0.173
    mu_pvdf_k = 0.154
    e_high, e_low, tau, alpha = 4, 54.2, 2.82e-6, 0.438
    e_cc_max = np.real(e_high + (e_low - e_high)/(1 + np.power(1j*2*np.pi*frequency*tau, 1 - alpha)))
    m_pvdf = 1.78e3*t_dielectric*L_dielectric*w_dielectric
    period = 1/2/frequency

    h2over4 = lambda h: np.square(h)/4
    h2 = lambda h: np.square(h)
    F32 = lambda h: np.exp(-h2over4(h))*np.sqrt(h)*((h2(h) + 1)*kv(0.25, h2over4(h)) - h2(h)*kv(0.75, h2over4(h)))/4/np.sqrt(np.pi)

    sigma_gw, C_gw = 2.76314774e-06, 3.05500533e+14
    F_es = lambda k, t_air, V: Fes_scalar/8*eps0*np.square(k)*(V**2)*(N*w_trace*w_pin)/np.sqrt(2*np.pi)/sigma_gw* \
                               quad(lambda y: np.divide(np.exp(-np.square(t_air - y)/2/np.square(sigma_gw)), np.square(t_dielectric + k*y)), -t_dielectric/k*0.83, 1e-3)[0]
    F_gw = lambda h: C_gw*np.power(sigma_gw, 1.5)*(N*w_trace*w_pin)*F32(h/sigma_gw)  # L_dielectric

    t_air_with_Fes = root_scalar(lambda t_air: F_gw(t_air) - F_es(k=e_cc_max, t_air=t_air, V=V_max) - m_pvdf*g - Fpreload, bracket=(1e-6, 10e-6)).root
    t_air_without_Fes = root_scalar(lambda t_air: F_gw(t_air) - m_pvdf*g - Fpreload, x0=3e-6, x1=4e-6).root

    string_to_write = "V = {}, f = {}Hz, Fpreload = {:0.3f} --> t_air, with Fes={:0.3e}, without Fes={:0.3e} --> Fpred = {:0.3f}, Ppred = {:0.3f} kPa"
    print(string_to_write.format(V, frequency, Fpreload, t_air_with_Fes, t_air_without_Fes,
                                 mu_pvdf_s*(F_gw(t_air_with_Fes) - F_gw(t_air_without_Fes)),
                                 mu_pvdf_s*(F_gw(t_air_with_Fes) - F_gw(t_air_without_Fes))/(w_pin*L_dielectric)/1e3))

    return (mu_pvdf_s*F_gw(t_air_with_Fes) - mu_pvdf_k*F_gw(t_air_without_Fes))/(w_pin*L_dielectric)/1e3


def fit_Fes_scalar(freq, Fes_scalar):
    return [calculate_shear_force_pred(V_max=V, Fpreload=F_preload, frequency=f, Fes_scalar=Fes_scalar) for f in freq]


all_freq = sorted(all_load_cell_max_force.keys())
ax2_lines = []
markers = ['o', '^', 'd', 'x']
colors = ["#1964B0", "#F1932D", "#4DB264", "#DB060B"]
all_x2, all_y2, all_yerr2 = [], [], []
best_Fes_scalars = []
best_Rleakages = []
F_preload = np.mean([np.mean([np.mean(all_load_cell_preload_forces[f][V]) for V in all_V]) for f in all_freq])
for i, V in enumerate(all_V):
    y = []
    yerr = []
    ypred = []
    for f in sorted(all_load_cell_max_force.keys()):
        y.append(np.mean(all_load_cell_max_force[f][V]))
        yerr.append(np.std(all_load_cell_max_force[f][V]))

    popt, pcov = curve_fit(fit_Fes_scalar, all_freq, y, p0=(1,), bounds=(0.1, 10))
    best_Fes_scalars.append(popt[0])
    print("Optimal Fes_scalar for V =", V, "=", popt[0])
    f_sim = np.power(10, np.linspace(-1, 4, 10))
    for f in f_sim:
        ypred.append(calculate_shear_force_pred(V_max=V, Fpreload=F_preload, frequency=f, Fes_scalar=popt[0]))
    line = ax2.errorbar(all_freq, y=y, yerr=yerr, fmt='-' + markers[i], lw=1, capsize=6, capthick=1, ms=6, label="{} V".format(int(V)),
                        c=colors[i])
    model_line, = ax2.plot(f_sim, ypred, ls='--', c=colors[i], label="Model")
    ax2_lines.append(line)

    all_x2.append(all_freq)
    all_y2.append(y)
    all_yerr2.append(yerr)
ax2_lines.append(model_line)
print("Best Fes scalars", best_Fes_scalars, np.mean(best_Fes_scalars), np.std(best_Fes_scalars))

# ax2.grid(True, which='minor', c='0.9')
ax2.grid(True)
ax2.set_xlabel("Drive Frequency (Hz)")
ax2.set_ylabel("Shear Force Capacity (kPa)")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xticks([0.1, 1, 10, 100, 1000, 10000], ["DC", r"$10^0$", r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$"])
ax2.set_yticks([1, 2, 5, 10], [1, 2, 5, 10])

ax2.legend([ax2_lines[3], ax2_lines[2], ax2_lines[1], ax2_lines[0], ax2_lines[4]], ["300 V", "250 V", "200 V", "150 V", "Model"],
           loc='lower left', ncol=3, fontsize=14)
fig2.text(0.00, 1, "(a)", transform=fig2.transFigure, horizontalalignment='left', verticalalignment='top')

# np.save(save_folder + timestamp + "_plotdata",  [all_V, all_freq, all_x1, all_y1, all_yerr1, all_x2, all_y2, all_yerr2,
#                                                  all_load_cell_disengage_forces, all_load_cell_preload_forces])
plt.savefig(save_folder + timestamp + ".png", dpi=300)
# # plt.savefig("figures/" + timestamp + ".svg")
plt.savefig(save_folder + timestamp + ".pdf")

plt.show()
