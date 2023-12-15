from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

RECONSTRUCTIONS = ["splinempe", "millipede"]
SPLINEMPE_INDEX = 0
MILLIPEDE_INDEX = 1
SPLINEMPE_FILENAME = "gcn_notices_gold_bronze.txt"
# IC190819A: Updated GCN Notice never reported
# IC191231A: in GCN Circular reported the day after
# IC200227A: Updated GCN Notice never reported
# IC210503A: Never reconstructed with SplineMPE because detector
# was in a test run configuration. Missing in SplineMPE catalog
# IC210510A: Updated GCN Notice never reported
# IC211117A: Updated GCN Notice never reported
# IC221223A: Updated GCN Notice never reported
# IC230120A: Updated GCN Notice never reported
SPLINEMPE_EXCEPTIONS = [
    "IC191231A",
    "IC190819A",
    "IC200227A",
    "IC210510A",
    "IC211117A",
    "IC221223A",
    "IC230220A",
]
SPLINEMPE_MISSING = "IC210510A"
# IC200120A: likely background
# IC230823A: likely background
SPLINEMPE_BACKGROUND = ["IC200120A", "IC230823A"]
SPLINEMPE_GCN_START = "<tr align=left>\n"
SPLINEMPE_INDEX_START = 65
SPLINEMPE_COMMENT_START = "<!--\n"
MILLIPEDE_FILENAME = "IC_Alerts_Table.csv"
MILLIPEDE_IC_NAME = "IC_NAME"
MILLIPEDE_RA = "RA"
MILLIPEDE_RA_PLUS = "RA_ERR_P"
MILLIPEDE_RA_MINUS = "RA_ERR_M"
MILLIPEDE_DEC = "DEC"
MILLIPEDE_DEC_PLUS = "DEC_ERR_P"
MILLIPEDE_DEC_MINUS = "DEC_ERR_M"
FIGSIZE = (5, 4)
RANGE_BINS = (-0.5, 3)
NBINS = 10
ALPHA = 0.8
HISTTYPE = "stepfilled"
LINEWIDTH = 2
EDGECOLOR = "darkblue"
LABEL_HIST = r"$\dfrac{\mathrm{Millipede}}{\mathrm{SplineMPE}}$"
AVERAGE_COLOR = "black"
AVERAGE_LINEWIDTH = 3
AVERAGE_LINESTYLE = "--"
LEGEND_FONTSIZE = "x-large"
TICKS_FONTSIZE = "x-large"
XLABEL = "Ratio"
YLABEL = "Counts"
AXIS_FONTSIZE = "xx-large"
TITLE = "Uncertainty areas"
TITLE_FONTSIZE = "xx-large"
PLOTNAME = "areas_ratio"

print("Definition of paths...")

# Definition of paths
cwd = Path(os.getcwd())
data_path = cwd / "../data"
figures_path = cwd / "../figures"

print(f"Retrieving the alerts reconstructed with {RECONSTRUCTIONS[SPLINEMPE_INDEX]}...")

NAMEs_spli = np.array([])
spline_areas = np.array([])
rev0 = False
rev1 = False
has_rev1 = False
is_data = False
rev1_names = np.array([])
splinempe_f = open(data_path / SPLINEMPE_FILENAME)
for index, line in enumerate(splinempe_f):
    if line == SPLINEMPE_GCN_START and index > SPLINEMPE_INDEX_START:
        notice_line_index = 0
        if index < 100:
            is_data = True
    if line == SPLINEMPE_COMMENT_START:
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
                    name in rev1_names or name in SPLINEMPE_EXCEPTIONS
                ) and not name in SPLINEMPE_BACKGROUND:
                    NAMEs_spli = np.append(NAMEs_spli, name)
                    has_rev1 = True
                else:
                    has_rev1 = False
        if rev0 and is_data and has_rev1:
            if notice_line_index == 9:
                err_90 = float(line.split(">")[1].split("<")[0]) / 60
                area90 = np.pi * np.deg2rad(err_90) ** 2
                spline_areas = np.append(spline_areas, area90)

print(f"Retrieving the alerts reconstructed with {RECONSTRUCTIONS[MILLIPEDE_INDEX]}...")

alerts_df = pd.read_csv(data_path / MILLIPEDE_FILENAME)
missing_index = alerts_df[alerts_df[MILLIPEDE_IC_NAME] == SPLINEMPE_MISSING].index[0]
alerts_df = alerts_df.drop(missing_index)
RAs = alerts_df[MILLIPEDE_RA].to_numpy()
DECs = alerts_df[MILLIPEDE_DEC].to_numpy()
RAs_ERR_PLUS = alerts_df[MILLIPEDE_RA_PLUS].to_numpy()
DECs_ERR_PLUS = alerts_df[MILLIPEDE_DEC_PLUS].to_numpy()
RAs_ERR_MINUS = alerts_df[MILLIPEDE_RA_MINUS].to_numpy()
DECs_ERR_MINUS = alerts_df[MILLIPEDE_DEC_MINUS].to_numpy()
NAMEs = alerts_df[MILLIPEDE_IC_NAME].to_numpy()


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

fig, ax = plt.subplots(figsize=FIGSIZE)
logbins = np.logspace(RANGE_BINS[0], RANGE_BINS[1], NBINS)
plt.hist(
    ratios,
    bins=logbins,
    alpha=ALPHA,
    histtype=HISTTYPE,
    linewidth=LINEWIDTH,
    edgecolor=EDGECOLOR,
    label=LABEL_HIST,
)
plt.axvline(
    np.mean(ratios),
    color=AVERAGE_COLOR,
    linewidth=AVERAGE_LINEWIDTH,
    linestyle=AVERAGE_LINESTYLE,
    label=f"Average: {str(int(np.mean(ratios)))}",
)
plt.legend(fontsize=LEGEND_FONTSIZE)
plt.xlim(0.24, 2000)
# plt.ylim(0,25)
ax.set_xticklabels(ax.get_xticks(), size=TICKS_FONTSIZE)
ax.set_yticklabels(ax.get_yticks(), size=TICKS_FONTSIZE)
plt.xscale("log")
plt.xlabel(XLABEL, size=AXIS_FONTSIZE)
plt.ylabel(YLABEL, size=AXIS_FONTSIZE)
plt.title(TITLE, size=TITLE_FONTSIZE)
plt.savefig(figures_path / PLOTNAME, bbox_inches="tight")
