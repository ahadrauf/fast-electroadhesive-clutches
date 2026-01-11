from sklearn.metrics import r2_score
from sim_engagement_time import *

now = datetime.now()
name_clarifier = "_engagement_time_vs_frequency_plotandsim"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
save_folder = "C:/Users/ahadrauf/Desktop/Research/latex/electroadhesive_dynamics_paper/figures_test/"

plt.rcParams["font.family"] = "Arial"
plt.rc('font', size=16)  # controls default text size
plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=16)  # fontsize of the y tick labels
plt.rc('legend', fontsize=16)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

file_loc = "../data/"
file_name = "20241208_14_37_18_timing_vs_frequency_plotandsim_plotdata"

data = np.load(file_loc + file_name + ".npy", allow_pickle=True)
all_V, all_x1, all_y1, all_yerr1, all_x2, all_y2, all_yerr2, all_load_cell_engage_times, all_load_cell_preload_forces, all_load_cell_initial_forces, all_load_cell_max_forces, all_load_cell_disengage_times, all_load_cell_disengage_forces, all_load_cell_disengage_times_10pct = data

all_freq = sorted(all_load_cell_preload_forces[150].keys())

num_readings = sum([sum([len(all_load_cell_engage_times[V][freq]) for freq in all_freq]) for V in all_V])
num_zeros = sum([sum([len(all_load_cell_engage_times[V][freq][np.where(all_load_cell_engage_times[V][freq] == 0)]) for freq in all_freq]) for V in all_V])
num_over_30 = sum([sum([len(all_load_cell_engage_times[V][freq][np.where(all_load_cell_engage_times[V][freq] > 30)]) for freq in all_freq]) for V in all_V])
for V in all_V:
    for freq in all_freq:
        for val in all_load_cell_engage_times[V][freq]:
            if val > 30:
                print("V = {}, freq = {} Hz, val = {}".format(V, freq, val))
print("Num Readings:", num_readings, "Num Zeros:", num_zeros, "Num Over 30:", num_over_30)

# Plot data
colors = ["#1964B0", "#F1932D", "#4DB264", "#DB060B"]  # "#F7F057" = yellow
fig, axs_all = plt.subplots(2, 2, layout='constrained', sharex='all', sharey='all',
                            figsize=(6, 4))  # , figsize=(10, 8)
axs = [axs_all[0][0], axs_all[0][1], axs_all[1][0], axs_all[1][1]]
for i, V in enumerate(all_V):
    x = all_x1[i]
    y = all_y1[i]
    yerr = all_yerr1[i]
    errorbar = axs[i].errorbar(x, y, yerr=yerr,  # xerr=[np.std(fi[i]) for fi in forces_by_bucket],
                               capsize=5, ecolor='k', elinewidth=2, capthick=2, c='k')  # alpha=0.4
    scatter = axs[i].scatter(x, y, c='k', s=50)
    print("V =", V, "Force:", y)

    frequency_range = np.power(10, np.linspace(-1, 4, 200))

    freq_sim = []
    tengage_sim = []
    all_preloads = np.hstack([all_load_cell_preload_forces[V][f] for f in all_freq])
    print("Preloads for V =", V, np.mean(all_preloads))
    Fpreload = np.mean(all_preloads)
    for frequency in frequency_range:
        L_dielectric = 55.5e-3
        w_pin = 2e-3
        depth_pin = 2e-3
        k_dielectric = 54.2
        t_dielectric = 24e-6
        T, Y, dY, X, dX, events, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic = sim_engagement_time(L_dielectric=L_dielectric,
                                                                                                               w_pin=w_pin, depth_pin=depth_pin,
                                                                                                               V_max=V, t_dielectric=t_dielectric,
                                                                                                               period=1 / 2 / frequency, Fext=Fpreload,
                                                                                                               k_dielectric=k_dielectric)
        if type(events) == np.float64:
            freq_sim.append(frequency)
            tengage_sim.append(events)
    axs[i].semilogx(freq_sim, tengage_sim, ls='--', c='tab:red', lw=2)

    bbox = dict(alpha=0.8, boxstyle='round', fc='white', ec='0.8')
    axs[i].text(0.05, 0.925, "{} V".format(int(V)), bbox=bbox, transform=axs[i].transAxes, horizontalalignment='left', verticalalignment='top')
    axs[i].grid(True)
    axs[i].set_xscale('log')
    axs[i].tick_params(labelleft=True)
    axs[i].set_xticks([0.1, 1, 10, 100, 1000, 10000], ["DC", r"$10^0$", r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$"], fontsize=16)
    axs[i].set_yticks([0, 10, 20])
    axs[i].xaxis.set_tick_params(labelbottom=True)

    tengage_r2 = []
    for frequency in x:
        L_dielectric = 55.5e-3
        w_pin = 2e-3
        depth_pin = 2e-3
        k_dielectric = 54.2
        t_dielectric = 24e-6
        T, Y, dY, X, dX, events, Fk, Fb, Fes, Fconstant, Flinear_static, Flinear_kinetic = sim_engagement_time(L_dielectric=L_dielectric,
                                                                                                               w_pin=w_pin, depth_pin=depth_pin,
                                                                                                               V_max=V, t_dielectric=t_dielectric,
                                                                                                               period=1 / 2 / frequency, Fext=Fpreload,
                                                                                                               k_dielectric=k_dielectric)
        tengage_r2.append(events)
    print("R2 score for V =", V, r2_score(y, tengage_r2))


fig.text(0.00, .98, "(c)", transform=fig.transFigure, horizontalalignment='left', verticalalignment='top')
fig.supxlabel("Drive Frequency (Hz)", fontsize=16)
fig.supylabel("Engagement Time (µs)", fontsize=16)

# plt.savefig(save_folder + timestamp + ".png", dpi=300)
# # # plt.savefig("figures/" + timestamp + ".svg")
# plt.savefig(save_folder + timestamp + ".pdf")

plt.show()
