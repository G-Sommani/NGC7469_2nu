from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import config as cfg

print("Definition of paths...")

# Definition of paths
cwd = Path(os.getcwd())
data_path = cwd / "../data"
figures_path = cwd / "../figures"

print(
    f"Retrieving the alerts reconstructed with {cfg.ALLOWED_RECONSTRUCTIONS[cfg.SPLINEMPE_INDEX]}..."
)

NAMEs_spli = np.array([])
spline_areas = np.array([])
rev0 = False
rev1 = False
has_rev1 = False
is_data = False
rev1_names = np.array([])
splinempe_f = open(data_path / cfg.SPLINEMPE_FILENAME)
for index, line in enumerate(splinempe_f):
    if line == cfg.SPLINEMPE_GCN_START and index > cfg.SPLINEMPE_INDEX_START:
        notice_line_index = 0
        if index < 100:
            is_data = True
    if line == cfg.SPLINEMPE_COMMENT_START:
        is_data = False
    if is_data:
        if notice_line_index == 2:
            rev_number = line.split(">")[1].split("<")[0]
            if rev_number == "1" or rev_number == "2":
                rev0 = False
                rev1 = True
            if rev_number == "0":
                rev0 = True
                rev1 = False
        notice_line_index += 1
        if notice_line_index == 4:
            date = line.split(">")[1].split("<")[0]
            year = date.split("/")[0]
            month = date.split("/")[1]
            day = date.split("/")[2]
            name = f"IC{year}{month}{day}A"
            if rev1:
                if len(rev1_names) > 0:
                    if name in rev1_names:
                        if name == f"IC{year}{month}{day}B":
                            continue
                        rev1_names[
                            np.where(rev1_names == name)
                        ] = f"IC{year}{month}{day}B"
                rev1_names = np.append(rev1_names, name)
            elif rev0:
                if len(NAMEs_spli) > 0:
                    if name in NAMEs_spli:
                        NAMEs_spli[
                            np.where(NAMEs_spli == name)
                        ] = f"IC{year}{month}{day}B"
                if (
                    name in rev1_names or name in cfg.SPLINEMPE_EXCEPTIONS
                ) and not name in cfg.SPLINEMPE_BACKGROUND:
                    NAMEs_spli = np.append(NAMEs_spli, name)
                    has_rev1 = True
                else:
                    has_rev1 = False
        if rev0 and is_data and has_rev1:
            if notice_line_index == 9:
                err_90 = float(line.split(">")[1].split("<")[0]) / 60
                area90 = np.pi * np.deg2rad(err_90) ** 2
                spline_areas = np.append(spline_areas, area90)

print(
    f"Retrieving the alerts reconstructed with {cfg.ALLOWED_RECONSTRUCTIONS[cfg.MILLIPEDE_INDEX]}..."
)

alerts_df = pd.read_csv(data_path / cfg.MILLIPEDE_FILENAME)
missing_index = alerts_df[
    alerts_df[cfg.MILLIPEDE_IC_NAME] == cfg.SPLINEMPE_MISSING
].index[0]
alerts_df = alerts_df.drop(missing_index)
RAs = alerts_df[cfg.MILLIPEDE_RA].to_numpy()
DECs = alerts_df[cfg.MILLIPEDE_DEC].to_numpy()
RAs_ERR_PLUS = alerts_df[cfg.MILLIPEDE_RA_PLUS].to_numpy()
DECs_ERR_PLUS = alerts_df[cfg.MILLIPEDE_DEC_PLUS].to_numpy()
RAs_ERR_MINUS = alerts_df[cfg.MILLIPEDE_RA_MINUS].to_numpy()
DECs_ERR_MINUS = alerts_df[cfg.MILLIPEDE_DEC_MINUS].to_numpy()
NAMEs = alerts_df[cfg.MILLIPEDE_IC_NAME].to_numpy()


def millipede_area(index):
    """
    Given the index of the alert, returns the angular area
    """
    phi1_deg = RAs[index] - RAs_ERR_MINUS[index]
    phi2_deg = RAs[index] + RAs_ERR_PLUS[index]
    delta1_deg = DECs[index] - DECs_ERR_MINUS[index]
    delta2_deg = DECs[index] + DECs_ERR_PLUS[index]
    phi1_rad = np.deg2rad(phi1_deg)
    phi2_rad = np.deg2rad(phi2_deg)
    delta1_rad = np.deg2rad(delta1_deg)
    delta2_rad = np.deg2rad(delta2_deg)
    A = (phi2_rad - phi1_rad) * (np.sin(delta2_rad) - np.sin(delta1_rad))
    return A


ratios = np.array([])
for i in range(len(alerts_df)):
    area = millipede_area(i)
    ratio = area / spline_areas[-i - 1]
    ratios = np.append(ratios, ratio)

print("Plotting...")

fig, ax = plt.subplots(figsize=cfg.FIGSIZE_AREAS)
logbins = np.logspace(cfg.RANGE_BINS_AREAS[0], cfg.RANGE_BINS_AREAS[1], cfg.NBINS_AREAS)
plt.hist(
    ratios,
    bins=logbins,
    alpha=cfg.ALPHA_AREAS,
    histtype=cfg.HISTTYPE_AREAS,
    linewidth=cfg.LINEWIDTH_AREAS,
    edgecolor=cfg.EDGECOLOR_AREAS,
    label=cfg.LABEL_HIST_AREAS,
)
plt.axvline(
    np.mean(ratios),
    color=cfg.AVERAGE_COLOR_AREAS,
    linewidth=cfg.AVERAGE_LINEWIDTH_AREAS,
    linestyle=cfg.AVERAGE_LINESTYLE_AREAS,
    label=f"Average: {str(int(np.mean(ratios)))}",
)
plt.legend(fontsize=cfg.LEGEND_FONTSIZE_AREAS)
plt.xlim(0.24, 2000)
# plt.ylim(0,25)
ax.set_xticklabels(ax.get_xticks(), size=cfg.TICKS_FONTSIZE_AREAS)
ax.set_yticklabels(ax.get_yticks(), size=cfg.TICKS_FONTSIZE_AREAS)
plt.xscale("log")
plt.xlabel(cfg.XLABEL_AREAS, size=cfg.AXIS_FONTSIZE_AREAS)
plt.ylabel(cfg.YLABEL_AREAS, size=cfg.AXIS_FONTSIZE_AREAS)
plt.title(cfg.TITLE_AREAS, size=cfg.TITLE_FONTSIZE_AREAS)
plt.savefig(figures_path / cfg.PLOTNAME_AREAS, bbox_inches="tight")
plotname_pdf = cfg.PLOTNAME_AREAS + ".pdf"
plt.savefig(figures_path / plotname_pdf, bbox_inches="tight")
