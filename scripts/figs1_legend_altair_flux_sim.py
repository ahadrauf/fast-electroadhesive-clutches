# This is a simple script to generate the legends I put on Fig. S1 (I then added these legends to the figures in PowerPoint)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from datetime import datetime

now = datetime.now()
name_clarifier = "_timing_lit_review"
timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
plt.rcParams["font.family"] = "Arial"
plt.rc('font', size=24)  # controls default text size
plt.rc('axes', labelsize=24)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=24)  # fontsize of the x tick labels
plt.rc('ytick', labelsize=24)  # fontsize of the y tick labels
plt.rc('legend', fontsize=24)  # fontsize of the legend
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

colors_top = list(reversed(["#ffffb0", "#fff700", "#ffdb00", "#ffb400", "#ff8000", "#ff3316",
                            "#f00050", "#cb008e", "#9600c5", "#4900e7", "#0000ac", "#000070"]))
top_max = 27.374
colors_bottom = list(reversed(["#efefef", "#dfdfdf", "#cfcfcf", "#bfbfbf", "#afafaf", "#9f9f9f",
                               "#8f8f8f", "#7f7f7f", "#6f6f6f", "#5f5f5f", "#4f4f4f", "#3f3f3f",
                               "#2f2f2f", "#1f1f1f", "#0f0f0f"]))
ticks_bottom = [0, 40, 80, 120, 160, 200, 240, 280]


def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

fig, ax = plt.subplots(1, 1, layout='constrained')
Z = top_max*np.random.rand(100, 100)
cmap = LinearSegmentedColormap.from_list(r"$E_d$ (V/m)", colors_top, N=100)
im = ax.imshow(Z, origin='upper', cmap=cmap)
im.set_clim([0, top_max])
cbar = fig.colorbar(im, ax=ax)  # , format=mticker.FuncFormatter(fmt))
# cbar = fig.colorbar(im, ax=ax, format="%4.1e")
cbar.set_label(r"Electric Field $E_{air}$ (V/μm)")

# fig, ax = plt.subplots(1, 1, layout='constrained')
# Z = 300*np.random.rand(100, 100)
# cmap = LinearSegmentedColormap.from_list("Potential (V)", colors_bottom, N=15)
# im = ax.imshow(Z, origin='upper', cmap=cmap)
# im.set_clim([0, 300])
# cbar = fig.colorbar(im, ax=ax, ticks=ticks_bottom)
# cbar.set_label("Electric Potential (V)")

plt.show()
