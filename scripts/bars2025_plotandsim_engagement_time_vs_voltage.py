import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit, root_scalar
from scipy.integrate import solve_ivp
from scipy.special import kv
from scipy.stats import linregress
from datetime import datetime
import csv
# from sim_engagement_time_20241121 import *
from sim_engagement_time import *

now = datetime.now()
name_clarifier = "_engagement_time_vs_voltage_plotandsim_bars2025"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
save_folder = "../figures/"

plt.rcParams["font.family"] = "Arial"
# plt.rcParams["font.family"] = "serif"
# plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rc('font', size=26)  # controls default text size
plt.rc('axes', labelsize=26)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=26)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=26)  # fontsize of the y tick labels
plt.rc('legend', fontsize=20)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
# plt.style.use('tableau-colorblind10')

file_loc = "../data/"
file_name = "20241105_11_23_31_timing_vs_voltage_plotdata"
file_name = "20241206_23_08_14_timing_vs_voltage_plotdata"

data = np.load(save_folder + file_name + ".npy", allow_pickle=True)
# all_V, x1, y1, yerr1, x2, y2, yerr2, all_load_cell_engage_times, all_load_cell_disengage_times, all_load_cell_disengage_forces, all_load_cell_preload_forces = data
all_V, x1, y1, yerr1, x2, y2, yerr2, all_load_cell_engage_times, all_load_cell_preload_forces, all_load_cell_initial_forces, all_load_cell_max_forces, all_load_cell_disengage_times, all_load_cell_disengage_forces, all_load_cell_disengage_times_10pct = data
print("Min Release Time", np.min(y2), "ms")
print("Average Preload Force", [np.mean(all_load_cell_preload_forces[V]) for V in all_V], "N")

slope, intercept, r, p, se = linregress(x1, y1)
print("Linear fit data:", slope, intercept, r*r)
print(x1)
print(y1)

# all_freq = sorted(all_load_cell_preload_forces.keys())
num_readings = sum([len(all_load_cell_engage_times[freq]) for freq in all_V])
num_zeros = sum([len(all_load_cell_engage_times[freq][np.where(all_load_cell_engage_times[freq] == 0)]) for freq in all_V])
num_over_30 = sum([len(all_load_cell_engage_times[freq][np.where(all_load_cell_engage_times[freq] > 30)]) for freq in all_V])
for freq in all_V:
    print("freq = {} V".format(freq), all_load_cell_engage_times[freq], np.mean(all_load_cell_engage_times[freq]))
    for val in all_load_cell_engage_times[freq]:
        if val > 30:
            print("freq = {} Hz, val = {}".format(freq, val))
print("Num Readings:", num_readings, "Num Zeros:", num_zeros, "Num Over 30:", num_over_30)

# Plot data
colors = ["#1964B0", "#F1932D", "#4DB264", "#DB060B"]  # "#F7F057" = yellow
# fig2, ax2 = plt.subplots(1, 1, layout='constrained', figsize=(6, 3.25))
fig2, ax2 = plt.subplots(1, 1, layout='constrained', figsize=(5.5, 4))
errorbar2 = ax2.errorbar(x1, y1, yerr=yerr1, capsize=5, ecolor='k', elinewidth=2, capthick=2, color='k', lw=2)
scatter2 = ax2.scatter(x1, y1, c='k', s=50)
ax2.grid(True)
ax2.set_xlabel("Voltage (V)")
ax2.set_ylabel("Engagement Time (µs)", y=0.4)
bbox = dict(alpha=0.8, boxstyle='round', fc='white', ec='0.8')
# ax2.text(0.96, 0.95, r"$w_s$" + " = 2 mm\nf = 1000 Hz", bbox=bbox, transform=ax2.transAxes, horizontalalignment='right', verticalalignment='top')

V_range = np.linspace(100, 300, 100)
V_sim = []
trelease_sim = []
all_preload_forces = np.hstack([all_load_cell_preload_forces[V] for V in all_V])
# print([all_load_cell_preload_forces[V] for V in all_V])
# print(all_preload_forces)
# Fpreload = np.mean(all_preload_forces)
Fpreload = np.max([np.mean(all_load_cell_preload_forces[V]) for V in all_V])
print("Average preload force:", Fpreload)
# for V in V_range:
#     L_dielectric = 55.5e-3
#     w_pin = 2e-3
#     depth_pin = 2e-3
#     k_dielectric = 54.2
#     t_dielectric = 24e-6
#     frequency = 1000
#     period = 1 / 2 / frequency
#     # Fpreload = np.mean(all_load_cell_preload_forces[V])
#     T, Y, dY, X, dX, events, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic = sim_engagement_time(L_dielectric=L_dielectric,
#                                                                                                            w_pin=w_pin, depth_pin=depth_pin,
#                                                                                                            V_max=V, t_dielectric=t_dielectric,
#                                                                                                            k_dielectric=54.2,
#                                                                                                            period=period, Fext=Fpreload)
#     if type(events) == np.float64:
#         V_sim.append(V)
#         trelease_sim.append(events)

V_sim = [100.0, 102.02020202020202, 104.04040404040404, 106.06060606060606, 108.08080808080808, 110.1010101010101, 112.12121212121212, 114.14141414141415, 116.16161616161617, 118.18181818181819, 120.20202020202021, 122.22222222222223, 124.24242424242425, 126.26262626262627, 128.2828282828283, 130.3030303030303, 132.32323232323233, 134.34343434343435, 136.36363636363637, 138.3838383838384, 140.40404040404042, 142.42424242424244, 144.44444444444446, 146.46464646464648, 148.4848484848485, 150.50505050505052, 152.52525252525254, 154.54545454545456, 156.56565656565658, 158.5858585858586, 160.60606060606062, 162.62626262626264, 164.64646464646466, 166.66666666666669, 168.6868686868687, 170.70707070707073, 172.72727272727275, 174.74747474747477, 176.7676767676768, 178.7878787878788, 180.80808080808083, 182.82828282828285, 184.84848484848487, 186.8686868686869, 188.8888888888889, 190.90909090909093, 192.92929292929296, 194.94949494949498, 196.96969696969697, 198.989898989899, 201.010101010101, 203.03030303030303, 205.05050505050505, 207.07070707070707, 209.0909090909091, 211.11111111111111, 213.13131313131314, 215.15151515151516, 217.17171717171718, 219.1919191919192, 221.21212121212122, 223.23232323232324, 225.25252525252526, 227.27272727272728, 229.2929292929293, 231.31313131313132, 233.33333333333334, 235.35353535353536, 237.37373737373738, 239.3939393939394, 241.41414141414143, 243.43434343434345, 245.45454545454547, 247.4747474747475, 249.4949494949495, 251.51515151515153, 253.53535353535355, 255.55555555555557, 257.5757575757576, 259.5959595959596, 261.61616161616166, 263.6363636363636, 265.6565656565657, 267.67676767676767, 269.69696969696975, 271.7171717171717, 273.7373737373738, 275.75757575757575, 277.7777777777778, 279.7979797979798, 281.81818181818187, 283.83838383838383, 285.8585858585859, 287.8787878787879, 289.89898989898995, 291.9191919191919, 293.93939393939394, 295.95959595959596, 297.979797979798, 300.0]
trelease_sim = [2.428292184038009, 2.4229213620775347, 2.417803311804739, 2.412916966772522, 2.4082434741576466, 2.4037661679499567, 2.3994699575664615, 2.395340760379661, 2.391365940568458, 2.3875339577408203, 2.3838342622166953, 2.380257192980466, 2.3767902460175154, 2.3734233807662983, 2.370154621332613, 2.366977039651497, 2.3638842445448076, 2.360870333738632, 2.3579298490965552, 2.3550653551223557, 2.3522666050512107, 2.34952757073655, 2.346843931115974, 2.3442119537005, 2.3416281651394044, 2.3390893369702037, 2.3365924569191767, 2.33413480398684, 2.331714021399444, 2.3293273243863504, 2.32697243452549, 2.32464710561349, 2.3223494603234025, 2.320077641239876, 2.3178299074983872, 2.4044026603521904, 2.3945292860039844, 2.3849471854847817, 2.375645989318084, 2.3666156979639537, 2.3578466826667714, 2.349329681201636, 2.341055792421951, 2.3330164682174583, 2.325203503055477, 2.317609641898395, 2.3102278104231826, 2.3030496206217865, 2.2960686794300242, 2.289277803708006, 2.2826705847958815, 2.27624089440041, 2.2699828515287317, 2.2638904138710654, 2.2579590322901724, 2.2521821485676385, 2.246555404217423, 2.2410740566232517, 2.2357335378207117, 2.2305298558646225, 2.225458981512766, 2.220516221772433, 2.2156976545911014, 2.258259671767488, 2.2562934329934112, 2.2543288659028713, 2.252365634075492, 2.250428092816294, 2.2484689320598745, 2.246510663892133, 2.2445530999686634, 2.2425960600557118, 2.240639462971695, 2.238683506506336, 2.2367278196460694, 2.2347722648358075, 2.2328168050872295, 2.2308614087704752, 2.22890099326246, 2.2269399337919547, 2.2249814775578423, 2.223023301702945, 2.2210653787035746, 2.2191076867921242, 2.217150210938714, 2.2151929408413857, 2.2132358724014094, 2.211279005935256, 2.209322346458165, 2.2073670889837196, 2.2054145932849853, 2.2034625082162487, 2.2015107550875213, 2.1995591785792454, 2.197607810472588, 2.1956566838982248, 2.1937058085986787, 2.191755209609565, 2.1898049529389265, 2.1878550808590074]
print(V_sim)
print(trelease_sim)
print(ax2.get_ylim())
ax2.set_ylim(-4.646277307275232, 22.201907321532964)
xlim = ax2.get_xlim()
ax2.fill_between(ax2.get_xlim(), y1=ax2.get_ylim()[0], y2=15, color='#FFC000', alpha=0.4)
ax2.set_xlim(xlim)

line_model, = ax2.plot(V_sim, trelease_sim, ls='--', c=colors[-1], lw=2)
ax2.legend([(scatter2, errorbar2), line_model], ["Data", "Model"])
# ax2.text(0.005, .96, "(a)", transform=fig2.transFigure, horizontalalignment='left', verticalalignment='top')

slope, intercept, r, p, se = linregress(V_sim, trelease_sim)
print("Linear fit model:", slope, intercept, r*r)

# plt.savefig(save_folder + timestamp + ".png", dpi=300)
# plt.savefig(save_folder + timestamp + ".svg")
# plt.savefig(save_folder + timestamp + ".pdf")

plt.show()
