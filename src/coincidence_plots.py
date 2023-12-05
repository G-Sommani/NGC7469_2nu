import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import os

RA_INDEX = 0
DE_INDEX = 1
PLUS_INDEX = 1
MINUS_INDEX = 0
NGC7469_COORDS = [345.8156, 8.873861]
CIRCULAR_23_COORDS = [345.85, 9.14]
CIRCULAR_23_ERRS = [[1.04, 0.90], [0.76, 0.81]]
CIRCULAR_22_COORDS = [346.11, 8.91]
CIRCULAR_22_ERRS = [[1.33, 1.26], [1.01, 0.95]]
CIRCULAR_PLOT_FILENAME = "GCNCirculars_ngc7469"
NOTICE_22_COORDS = [345.757, 8.856]
NOTICE_22_ERR = 0.66
NOTICE_23_COORDS = [345.825, 9.007]
NOTICE_23_ERR = 0.513
NOTICE_PLOT_FILENAME = "GCNNotices_ngc7469"

print("Definition of paths...")

# Definition of paths
cwd = Path(os.getcwd())
data_path = cwd / "../data"
figures_path = cwd / "../figures"


def build_vertices(coords, errs):
    ra = coords[RA_INDEX]
    de = coords[DE_INDEX]
    ra_vertices = [
        ra - errs[RA_INDEX][MINUS_INDEX],
        ra - errs[RA_INDEX][MINUS_INDEX],
        ra + errs[RA_INDEX][PLUS_INDEX],
        ra + errs[RA_INDEX][PLUS_INDEX],
        ra - errs[RA_INDEX][MINUS_INDEX],
    ]
    de_vertices = [
        de - errs[DE_INDEX][MINUS_INDEX],
        de + errs[DE_INDEX][PLUS_INDEX],
        de + errs[DE_INDEX][PLUS_INDEX],
        de - errs[DE_INDEX][MINUS_INDEX],
        de - errs[DE_INDEX][MINUS_INDEX],
    ]
    return ra_vertices, de_vertices


def uncertainty_circle(coords, err):
    angles = np.linspace(0, 2 * np.pi, 100)
    ras = coords[RA_INDEX] + err * np.cos(angles)
    des = coords[DE_INDEX] + err * np.sin(angles)
    return ras, des


print("Plotting GCN Circulars...")

fig, ax = plt.subplots(figsize=(5, 5))

ra_vertices_22, de_vertices_22 = build_vertices(CIRCULAR_22_COORDS, CIRCULAR_22_ERRS)
ra_vertices_23, de_vertices_23 = build_vertices(CIRCULAR_23_COORDS, CIRCULAR_23_ERRS)

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
    CIRCULAR_22_COORDS[RA_INDEX],
    CIRCULAR_22_COORDS[DE_INDEX],
    marker="o",
    color="tab:blue",
)
plt.plot(
    CIRCULAR_23_COORDS[RA_INDEX],
    CIRCULAR_23_COORDS[DE_INDEX],
    color="tab:orange",
    marker="o",
)
plt.plot(
    NGC7469_COORDS[RA_INDEX],
    NGC7469_COORDS[DE_INDEX],
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

plt.savefig(figures_path / CIRCULAR_PLOT_FILENAME, bbox_inches="tight")
plt.close()

print("Plotting first GCN Notices...")

fig, ax = plt.subplots(figsize=(5, 5))

ras_err_22, des_err_22 = uncertainty_circle(NOTICE_22_COORDS, NOTICE_22_ERR)
ras_err_23, des_err_23 = uncertainty_circle(NOTICE_23_COORDS, NOTICE_23_ERR)

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
    NOTICE_22_COORDS[RA_INDEX], NOTICE_22_COORDS[DE_INDEX], marker="o", color="tab:blue"
)
plt.plot(
    NOTICE_23_COORDS[RA_INDEX],
    NOTICE_23_COORDS[DE_INDEX],
    color="tab:orange",
    marker="o",
)
plt.plot(
    NGC7469_COORDS[RA_INDEX],
    NGC7469_COORDS[DE_INDEX],
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

plt.savefig(figures_path / NOTICE_PLOT_FILENAME, bbox_inches="tight")
plt.close()
