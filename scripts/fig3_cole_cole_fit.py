import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from datetime import datetime

now = datetime.now()
name_clarifier = "_dielectric_vs_frequency_fit"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
plt.rcParams["font.family"] = "Arial"
plt.rc('font', size=20)  # controls default text size
plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=16)  # fontsize of the y tick labels
plt.rc('legend', fontsize=16)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

with open('../data/Wang 2009 - Dielectric Constant - v2.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    freqs_wang = []
    e_wang = []
    for row in reader:
        freqs_wang.append(float(row[0]))
        e_wang.append(float(row[1]))

freqs_bobnar = np.array([20, 100, 1e3, 10e3, 100e3, 1e6])
e_bobnar = np.array([40.814, 39.666, 37.683, 34.656, 25.678, 11.378])

freqs_polyk1 = np.array([100, 1e3, 10e3, 100e3, 1e6])
e_polyk1 = np.array([61.22, 57.018, 51.291, 36.4, 15.273])

freqs_polyk2 = np.array([1e3, 10e3, 100e3, 1e6])
e_polyk2 = np.array([53.196, 44.533, 29.662, 12.980])

freqs_kim = np.array([1e3, 10e3, 100e3])
e_kim = np.array([46.128, 41.815, 35.039])

with open('../data/Tang 2012 - Dielectric Constant.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    freqs_tang = []
    e_tang = []
    for row in reader:
        freqs_tang.append(float(row[0]))
        e_tang.append(float(row[1]))

freqs = np.concatenate((np.tile(freqs_wang, 5), np.tile(freqs_tang, 8), np.tile(freqs_polyk2, 101)))
e = np.concatenate((np.tile(e_wang, 5), np.tile(e_tang, 8), np.tile(e_polyk2, 101)))

# freqs_new = np.power(10, np.linspace(np.log10(np.min(freqs_wang)), np.log10(np.max(freqs_wang)), 30))
freqs_new = np.power(10, np.linspace(1, 11, 30))

with open('../data/cole_cole_inverse_laplace_pvdftrfecfe.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    t_cc = []
    e_cc = []
    for row in reader:
        t_cc.append(float(row[0]))
        e_cc.append(float(row[4]))
print("Min t:", 1 / max(freqs_wang), "Max t:", 1 / min(freqs_wang))

e_high = 4


def cole_cole(f, e_low, tau, alpha):
    return np.real(e_high + (e_low - e_high) / (1 + np.power(1j * 2 * np.pi * f * tau, alpha)))


def cole_cole_losstan(f, e_low, tau, alpha):
    return np.imag(e_high + (e_low - e_high) / (1 + np.power(1j * 2 * np.pi * f * tau, alpha)))


def cole_cole_inverse_laplace(t):
    return np.interp(t, t_cc, e_cc)


def havriliak_negami(f, e_low, tau, alpha, beta):
    return np.real(e_high + (e_low - e_high) / np.power(1 + np.power(1j * 2 * np.pi * f * tau, alpha), beta))


popt, pcov = curve_fit(cole_cole, freqs, e, p0=(50, 4e-6, 0.5), bounds=(0, (60, 100e-6, 1)), x_scale=(50, 1e-6, 1))
# popt, pcov = curve_fit(havriliak_negami, freqs, e, p0=(50, 4e-6, 0.5, 0.5),
#                        bounds=(0, (60, 100e-6, 1, 1)),
#                        x_scale=[50, 1e-6, 1, 1])

fig, ax1 = plt.subplots(1, 1, layout='constrained', figsize=(6, 4))
ax1.semilogx(1 / np.array(freqs_wang) / 2 / np.pi, e_wang, ls='--', c='tab:blue')
ax1.semilogx(1 / np.array(freqs_tang) / 2 / np.pi, e_tang, ls=':', c='tab:green')
ax1.semilogx(1 / freqs_polyk2 / 2 / np.pi, e_polyk2, ls='-.', c='tab:orange', zorder=4)
ax1.semilogx(1 / freqs_new, cole_cole_inverse_laplace(1 / freqs_new), c='k', lw=3)
ax1.set_xlabel("Time (s)")
ax1.set_ylabel(r"Relative Permittivity Re($\kappa$)")
ax1.legend(["Wang 2009", "Tang 2012", "Vendor Datasheet", "Cole-Cole Fit"])
ax1.set_xticks([1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1])
ax1.grid(True)

# plt.savefig("../figures/" + timestamp + ".png", dpi=300)
# # plt.savefig("figures/" + timestamp + ".svg")
# plt.savefig("../figures/" + timestamp + ".pdf")
plt.show()
