import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

now = datetime.now()
name_clarifier = "_timing_lit_review"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
plt.rcParams["font.family"] = "Arial"
plt.rc('font', size=16)  # controls default text size
plt.rc('axes', labelsize=16)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=16)  # fontsize of the y tick labels
plt.rc('legend', fontsize=16)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42


class PrevWork:
    def __init__(self, name, shear_forces=(), voltages=(), areas_cm2=(), shear_force_densities=(),
                 engagement_time=None, release_time=None, substrate=None, interdigitated_electrodes=True,
                 arrow_dx=10, arrow_dy=10, annotation_va="center", annotation_ha="left",
                 arrow_relpos=(0.5, 0.5)):
        self.name = name
        self.shear_forces = shear_forces  # N
        self.voltages = voltages
        self.areas_cm2 = areas_cm2
        self.shear_force_densities = shear_force_densities if shear_force_densities else \
            [shear_force / area * 10 for shear_force, area in zip(self.shear_forces, self.areas_cm2)]
        self.engagement_time = engagement_time
        self.release_time = release_time
        self.interdigitated_electrodes = interdigitated_electrodes
        self.substrate = substrate
        self.arrow_dx = arrow_dx
        self.arrow_dy = arrow_dy
        self.annotation_va = annotation_va
        self.annotation_ha = annotation_ha
        self.arrow_relpos = arrow_relpos

    def get_max_shear_force_density(self, return_idx=False):
        arr = self.shear_force_densities
        if return_idx:
            return np.max(arr), np.argmax(arr)
        else:
            return np.max(arr)

    def get_max_shear_force_density_per_volt(self, return_idx=False):  # kPa/V^2
        arr = [shear_force_density / volt for shear_force_density, volt in
               zip(self.shear_force_densities, self.voltages)]
        if return_idx:
            return np.max(arr), np.argmax(arr)
        else:
            return np.max(arr)

    def get_max_shear_force_density_per_volt2(self, return_idx=False):  # Pa/V^2
        arr = [shear_force_density / (volt**2) * 1e3 for shear_force_density, volt in
               zip(self.shear_force_densities, self.voltages)]
        if return_idx:
            return np.max(arr), np.argmax(arr)
        else:
            return np.max(arr)


lit = [PrevWork("This Work", shear_forces=(4.756,), voltages=(250,), areas_cm2=(1.11,), engagement_time=15e-6, release_time=874.6e-6, substrate="Brass",
                arrow_dx=24, arrow_dy=20, annotation_ha='center', annotation_va='bottom'),
       PrevWork("Zhang 2019", shear_forces=(0.946, 0.798, 0.5796, 0.3942, 0.290, 0.437, 0.374, 0.406, 0.0789, 0.8699, 0.8365, 0.7571, 0.77277, 0.6482, 0.4197),
                voltages=(355, 294, 250, 194, 160, 194, 160, 124, 84, 335, 294, 250, 194, 160, 124),
                areas_cm2=15 * [6.03 * 0.16], engagement_time=6.7e-3, substrate="Brass",
                arrow_dx=0, arrow_dy=-20, annotation_ha='center', annotation_va='top'),
       PrevWork("Shultz 2023", shear_forces=(0.0807, 0.933), voltages=(250 / np.sqrt(2), 250 / np.sqrt(2)), areas_cm2=(1, 1)),
       PrevWork("Kim 2023", shear_forces=(10.15, 7.02), voltages=(3500, 1500), areas_cm2=(4.776, 4.776)),
       PrevWork("Choi 2020", shear_force_densities=(23.9, 50, 10), voltages=(1000, 800, 100)),
       PrevWork("Cacucciolo 2019", shear_forces=(10,), voltages=(3500,), areas_cm2=(2,), release_time=150e-3, interdigitated_electrodes=True, substrate="PMMA Plastic",
                arrow_dx=-20, arrow_dy=0, annotation_ha='right', annotation_va='center', arrow_relpos=(0.5, 0.5)),  # previously, I said release time = 30e-3. Not sure why? 150ms based on Fig. 5a (many release curves have ~1.5 points, and force readings are sampled every 100 ms)),
       PrevWork("Wang 2014", shear_forces=(6.68,), voltages=(600 / np.sqrt(2),), areas_cm2=(27.4 * 17.3,)),
       PrevWork("Liu 2013", shear_forces=(37, 45.6), voltages=(2000, 3000), areas_cm2=(30 * 24, 30 * 24)),
       PrevWork("Yamamoto 2007", shear_forces=(18.85,), voltages=(1400 / np.sqrt(2),), areas_cm2=(7.5 * 13,)),
       PrevWork("Yehya 2016", shear_forces=(141.131, 74.1636, 33.55, 25.5, 9.516), voltages=(8000, 5000, 3000, 8000, 4000), areas_cm2=(400, 400, 400, 64, 64)),
       PrevWork("Choi 2019", shear_force_densities=(14.6, 10.8, 8.74), voltages=(1000, 1800, 2000)),
       PrevWork("Digumarti 2021", shear_force_densities=(0.534, 0.4222, 0.409), voltages=(1500, 1200, 1000)),
       PrevWork("Ruffatto 2014", shear_force_densities=(11.3, 45.9, 62, 36), voltages=(5000, 5000, 5000, 5000)),
       PrevWork("Prahlad 2008", shear_force_densities=(4.4, 2.1, 2.4, 4.1, 1.7, 0.8, 14), voltages=7 * [4000]),
       PrevWork("Cacucciolo 2022", shear_forces=(1.4059, 0.69037, 0.092998), voltages=(3000, 2000, 1000), areas_cm2=(2 * 3, 2 * 3, 2 * 3),
                engagement_time=300e-3, release_time=300e-3, interdigitated_electrodes=True, substrate="Fruit",
                arrow_dx=-20, arrow_dy=45, annotation_ha='center', annotation_va='bottom', arrow_relpos=(0.75, 0.5)),
       PrevWork("Hinchet 2020", engagement_time=5e-3, release_time=15e-3, interdigitated_electrodes=False,
                arrow_dx=10, arrow_dy=-30, annotation_ha='left', annotation_va='center', arrow_relpos=(0, 0.5)),  # release_time=8.7069e-3  (10% fall)
       PrevWork("Diller 2018", engagement_time=7.233e-3, release_time=23.6e-3, interdigitated_electrodes=False,
                arrow_dx=-15, arrow_dy=15, annotation_ha='right', annotation_va='center', arrow_relpos=(1, 0.5)),  # release_time=1.9531e-3  (10% fall)
       PrevWork("Wei 2023", engagement_time=50e-3, release_time=37e-3, interdigitated_electrodes=False,
                arrow_dx=20, arrow_dy=-20, annotation_ha='left', annotation_va='center', arrow_relpos=(0, 0.5)),  # unchanged release time based on condition
       PrevWork("Fitch 1957", engagement_time=150e-6, interdigitated_electrodes=False,
                arrow_dx=0, arrow_dy=-20, annotation_ha='center', annotation_va='top'),
       PrevWork("Gao 2019", release_time=36, interdigitated_electrodes=True, substrate="PVC",
                arrow_dx=-20, arrow_dy=0, annotation_ha='right', annotation_va='center'),
       PrevWork("Cao 2019", release_time=1.011e-1, interdigitated_electrodes=True, substrate="Cardboard",
                arrow_dx=-20, arrow_dy=-20, annotation_ha='right', annotation_va='top', arrow_relpos=(1, 1)),
       PrevWork("Yoder 2023", engagement_time=50e-3, release_time=70e-3, interdigitated_electrodes=False,
                arrow_dx=20, arrow_dy=0, annotation_ha='left', annotation_va='center', arrow_relpos=(0, 0.5)),  # release measured at 10% max force
       PrevWork("Chen 2022", release_time=5, interdigitated_electrodes=True, substrate="Glass",
                arrow_dx=-15, arrow_dy=15, annotation_ha='right', annotation_va='center', arrow_relpos=(1, 0.5)),
       PrevWork("Nakamura 2017", engagement_time=10.37e-3, release_time=161.6e-3, interdigitated_electrodes=False,
                arrow_dx=-12, arrow_dy=75, annotation_ha='center', annotation_va='bottom', arrow_relpos=(0.1, 0)),  # measured by exporting pdf into Illustrator (0.50432-0.49395), release to 0.10714+0.1*(0.85255-0.10714), 4.6562-4.4946
       PrevWork("Hinchet 2022", release_time=65e-3, interdigitated_electrodes=False,  # release_time=35.865e-3 (from sample plot)
                arrow_dx=-10, arrow_dy=-46, annotation_ha='right', annotation_va='top', arrow_relpos=(0.8, 0.5)),
       PrevWork("Jeon 1998", engagement_time=28.3e-3, release_time=13.69, interdigitated_electrodes=False, substrate="Glass",
                arrow_dx=5, arrow_dy=15, annotation_ha='left', annotation_va='center', arrow_relpos=(0.5, 0.5)),  # release time = 82.117e-3 if measured like with engagement time
       PrevWork("Li 2023", engagement_time=15.50e-3, release_time=43.47e-3, interdigitated_electrodes=False,
                arrow_dx=-15, arrow_dy=15, annotation_ha='right', annotation_va='center', arrow_relpos=(1, 0.5)),
       PrevWork("Mastrangelo 2023", engagement_time=419e-3, release_time=328e-3, interdigitated_electrodes=True, substrate="PLA Plastic",
                arrow_dx=17, arrow_dy=65, annotation_ha='center', annotation_va='center', arrow_relpos=(0.75, 0.5)),
       PrevWork("Thilakarathna 2024", engagement_time=270e-3, release_time=270e-3, interdigitated_electrodes=False,
                arrow_dx=-3, arrow_dy=33, annotation_ha='right', annotation_va='bottom', arrow_relpos=(0.75, 1)),
       PrevWork("Li 2021", engagement_time=25e-3, release_time=324e-3, interdigitated_electrodes=False,
                arrow_dx=0, arrow_dy=15, annotation_ha='center', annotation_va='bottom', arrow_relpos=(0.5, 0.5)),
       PrevWork("Krimsky 2024", engagement_time=24e-3, release_time=30e-3, interdigitated_electrodes=False,
                arrow_dx=12, arrow_dy=-30, annotation_ha='center', annotation_va='center', arrow_relpos=(0.34, 0.5)),
       PrevWork("Chen 1992", engagement_time=200e-3, release_time=10, interdigitated_electrodes=False,
                arrow_dx=20, arrow_dy=15, annotation_ha='left', annotation_va='center', arrow_relpos=(0, 0.5)),
       ]

fig, ax = plt.subplots(1, 1, figsize=(6, 5), layout='constrained')
all_markers = []
all_x = []
all_y = []
all_labels = []
all_annotation_params = []

for l in lit:
    if l.engagement_time is None and l.release_time is None:
        continue

    x1 = l.engagement_time
    y1 = l.release_time
    if x1 is None:
        x1 = 1e1
    if y1 is None:
        y1 = 1e2
    all_x.append(x1)
    all_y.append(y1)
    all_labels.append(l.name)
    all_markers.append('o')
    all_annotation_params.append((l.arrow_dx, l.arrow_dy, l.annotation_ha, l.annotation_va, l.arrow_relpos))

print("Number of papers:", len(lit))
print("Number of points to draw:", len(all_x))
xmin, xmax = 1e-5, 2e1
ymin, ymax = 5e-4, 2e2
offset_pct_x = lambda x, dpct: x * np.power(10, dpct * np.log10(xmax / xmin))
offset_pct_y = lambda y, dpct: y * np.power(10, dpct * np.log10(ymax / ymin))
ax.hlines(1e2, xmin=xmin, xmax=xmax, colors='r', alpha=0.4)
ax.vlines(10, ymin=ymin, ymax=ymax, colors='r', alpha=0.4)
scatter_plot = ax.scatter(all_x, all_y, marker='o', edgecolors='k', c='y', s=40, zorder=99)

count = 0
for x, y, label, annotation_params in zip(all_x, all_y, all_labels, all_annotation_params):
    dx, dy, ha, va, arrow_relpos = annotation_params
    color = 'k' if label == "This Work" else '0.3'
    ec = 'k' if label == 'This Work' else None
    visible = True if label == 'This Work' else False
    pad = 3 if label == 'This Work' else 0
    ax.annotate(label, xy=(x, y), xycoords='data', xytext=(dx, dy),
                textcoords='offset points', va=va, ha=ha, fontsize=12,
                arrowprops=dict(facecolor=color, arrowstyle='-|>', connectionstyle='arc3', relpos=arrow_relpos),
                bbox=dict(pad=pad, visible=visible, fill=False, fc=None, ec=ec),
                color=color)
    count += 1

ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.fill_between([xmin, 1e-3], ymin, 1e-3, color='tab:orange', alpha=0.3)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel("Engagement Time (s)")
ax.set_ylabel("Release Time (s)")
ax.set_xticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10],
              [r"$10^{-5}$", r"$10^{-4}$", r"$10^{-3}$", r"$10^{-2}$", r"$10^{-1}$", r"$10^{0}$", "DNM"])
ax.set_yticks([1e-3, 1e-2, 1e-1, 1, 10, 100],
              [r"$10^{-3}$", r"$10^{-2}$", r"$10^{-1}$", r"$10^{0}$", r"$10^{1}$", "DNM"])
plt.grid(True)

# plt.savefig("../figures/" + timestamp + ".png", dpi=300)
# # plt.savefig(save_folder + timestamp + ".svg")
# plt.savefig("../figures/" + timestamp + ".pdf")

# Automatically create lines to different scatter points: https://stackoverflow.com/a/74317297
plt.show()
