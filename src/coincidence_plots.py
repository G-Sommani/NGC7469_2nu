import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import os
import config as cfg

print("Definition of paths...")

# Definition of paths
cwd = Path(os.getcwd())
data_path = cwd / "../data"
figures_path = cwd / "../figures"


def build_vertices(coords, errs):
    ra_index = cfg.RA_INDEX
    de_index = cfg.DE_INDEX
    plus_index = cfg.PLUS_INDEX
    minus_index = cfg.MINUS_INDEX
    ra = coords[ra_index]
    de = coords[de_index]
    ra_vertices = [
        ra - errs[ra_index][minus_index],
        ra - errs[ra_index][minus_index],
        ra + errs[ra_index][plus_index],
        ra + errs[ra_index][plus_index],
        ra - errs[ra_index][minus_index],
    ]
    de_vertices = [
        de - errs[de_index][minus_index],
        de + errs[de_index][plus_index],
        de + errs[de_index][plus_index],
        de - errs[de_index][minus_index],
        de - errs[de_index][minus_index],
    ]
    return ra_vertices, de_vertices


def uncertainty_circle(coords, err):
    angles = np.linspace(0, 2 * np.pi, 100)
    ras = coords[cfg.RA_INDEX] + err * np.cos(angles)
    des = coords[cfg.DE_INDEX] + err * np.sin(angles)
    return ras, des


print("Plotting GCN Circulars...")

fig, ax = plt.subplots(figsize=(5, 5))

ra_vertices_22, de_vertices_22 = build_vertices(
    cfg.CIRCULAR_22_COORDS, cfg.CIRCULAR_22_ERRS
)
ra_vertices_23, de_vertices_23 = build_vertices(
    cfg.CIRCULAR_23_COORDS, cfg.CIRCULAR_23_ERRS
)

plt.plot(
    ra_vertices_22,
    de_vertices_22,
    color="tab:blue",
    linestyle="--",
    label="IC220424A\nGCN Circular",
    linewidth=3,
)
plt.plot(
    ra_vertices_23,
    de_vertices_23,
    color="tab:orange",
    linestyle="--",
    label="IC230416A\nGCN Circular",
    linewidth=3,
)
plt.plot(
    cfg.CIRCULAR_22_COORDS[cfg.RA_INDEX],
    cfg.CIRCULAR_22_COORDS[cfg.DE_INDEX],
    marker="o",
    color="tab:blue",
)
plt.plot(
    cfg.CIRCULAR_23_COORDS[cfg.RA_INDEX],
    cfg.CIRCULAR_23_COORDS[cfg.DE_INDEX],
    color="tab:orange",
    marker="o",
)
plt.plot(
    cfg.NGC7469_COORDS[cfg.RA_INDEX],
    cfg.NGC7469_COORDS[cfg.DE_INDEX],
    marker="*",
    markeredgecolor="black",
    color="yellow",
    markersize=20,
    label="NGC 7469",
    linestyle="",
    markeredgewidth=1.5,
)
plt.xlabel("Right Ascension [$^\circ$]", size="x-large")
plt.ylabel("Declination [$^\circ$]", size="x-large")
plt.xlim(344.5, 348)
plt.ylim(7.8, 11.5)
plt.title("GCN Circulars' 90% contours", size="xx-large")
plt.legend(fontsize="x-large")
ax.set_xticks([345, 346, 347, 348])
ax.set_xticklabels([345, 346, 347, 348], size="x-large")
ax.set_yticklabels(ax.get_yticks(), size="x-large")

print("Saving GCN Circulars plot...")

plt.savefig(figures_path / cfg.CIRCULAR_PLOT_FILENAME, bbox_inches="tight")
circular_plot_filename_pdf = cfg.CIRCULAR_PLOT_FILENAME + ".pdf"
plt.savefig(figures_path / circular_plot_filename_pdf, bbox_inches="tight")
plt.close()

print("Plotting first GCN Notices...")

fig, ax = plt.subplots(figsize=(5, 5))

ras_err_22, des_err_22 = uncertainty_circle(cfg.NOTICE_22_COORDS, cfg.NOTICE_22_ERR)
ras_err_23, des_err_23 = uncertainty_circle(cfg.NOTICE_23_COORDS, cfg.NOTICE_23_ERR)

plt.plot(
    ras_err_22,
    des_err_22,
    color="tab:blue",
    linestyle="--",
    label="IC220424A\nGCN Notice",
    linewidth=3,
)
plt.plot(
    ras_err_23,
    des_err_23,
    color="tab:orange",
    linestyle="--",
    label="IC230416A\nGCN Notice",
    linewidth=3,
)
plt.plot(
    cfg.NOTICE_22_COORDS[cfg.RA_INDEX],
    cfg.NOTICE_22_COORDS[cfg.DE_INDEX],
    marker="o",
    color="tab:blue",
)
plt.plot(
    cfg.NOTICE_23_COORDS[cfg.RA_INDEX],
    cfg.NOTICE_23_COORDS[cfg.DE_INDEX],
    color="tab:orange",
    marker="o",
)
plt.plot(
    cfg.NGC7469_COORDS[cfg.RA_INDEX],
    cfg.NGC7469_COORDS[cfg.DE_INDEX],
    marker="*",
    markeredgecolor="black",
    color="yellow",
    markersize=20,
    label="NGC 7469",
    linestyle="",
    markeredgewidth=1.5,
)
plt.xlabel("Right Ascension [$^\circ$]", size="x-large")
plt.ylabel("Declination [$^\circ$]", size="x-large")
plt.xlim(344.5, 348)
plt.ylim(7.8, 11.5)
plt.title("GCN Notices' 90% contours", size="xx-large")
plt.legend(fontsize="x-large")
ax.set_xticks([345, 346, 347, 348])
ax.set_xticklabels([345, 346, 347, 348], size="x-large")
ax.set_yticklabels(ax.get_yticks(), size="x-large")

print("Saving GCN Notices plot...")

plt.savefig(figures_path / cfg.NOTICE_PLOT_FILENAME, bbox_inches="tight")
notice_plot_filename_pdf = cfg.NOTICE_PLOT_FILENAME + ".pdf"
plt.savefig(figures_path / notice_plot_filename_pdf, bbox_inches="tight")
plt.close()
