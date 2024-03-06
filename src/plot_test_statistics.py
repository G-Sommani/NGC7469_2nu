import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path

SPLINEMPE = "splinempe"
MILLIPEDE = "millipede"
TURIN = "turin"
MILLIQUAS = "milliquas"
REDSHIFT = "redshift"
XRAY = "xray"
TS_FILENAME = "test_statistic"
FIGSIZE = (5, 3.5)
TURIN_LABEL = "Background distrubution,\nTurin catalog"
MILLIQUAS_LABEL = "Background distrubution,\nMilliquas catalog"
HISTTYPE = "step"
NBINS = 50
LINEWIDTH = 2
AXVCOLOR = "red"
AXVLINESTYLE = "--"
AXVLABEL = "IC220424A & IC230416A\nwith NGC 7469"
XLABEL = "Negative test statistic"
YLABEL = "N Datasets"
TITLE = "Goodness of Fit test"
XRAYTITLE = "Goodness of Fit test, X-ray weighting"
MILLIPEDETITLE = "Goodness of Fit test, Millipede"
MILLIPEDEXRAYTITLE = "Goodness of Fit test, Millipede, X-ray weighting"
FONTSIZE = "large"
AXISSCALE = "log"
FIGNAME = "TS_distr"
BBOX_INCHES = "tight"

print("Definition of paths...")

cwd = Path(os.getcwd())
data_results_path = cwd / "../data_results"
figures_path = cwd / "../figures"

print("Loading results...")

ts_spl_tur_z = np.load(data_results_path / f"{TS_FILENAME}_{SPLINEMPE}_{TURIN}_{REDSHIFT}.npy")
ts_spl_tur_z_result = np.load(data_results_path / f"{TS_FILENAME}_{SPLINEMPE}_{TURIN}_{REDSHIFT}_result.npy")

ts_spl_mil_z = np.load(data_results_path / f"{TS_FILENAME}_{SPLINEMPE}_{MILLIQUAS}_{REDSHIFT}.npy")
ts_spl_mil_z_result = np.load(data_results_path / f"{TS_FILENAME}_{SPLINEMPE}_{MILLIQUAS}_{REDSHIFT}_result.npy")

ts_spl_tur_xray = np.load(data_results_path / f"{TS_FILENAME}_{SPLINEMPE}_{TURIN}_{XRAY}.npy")
ts_spl_tur_xray_result = np.load(data_results_path / f"{TS_FILENAME}_{SPLINEMPE}_{TURIN}_{XRAY}_result.npy")

ts_mil_tur_z = np.load(data_results_path / f"{TS_FILENAME}_{MILLIPEDE}_{TURIN}_{REDSHIFT}.npy")
ts_mil_tur_z_result = np.load(data_results_path / f"{TS_FILENAME}_{MILLIPEDE}_{TURIN}_{REDSHIFT}_result.npy")

ts_mil_mil_z = np.load(data_results_path / f"{TS_FILENAME}_{MILLIPEDE}_{MILLIQUAS}_{REDSHIFT}.npy")
ts_mil_mil_z_result = np.load(data_results_path / f"{TS_FILENAME}_{MILLIPEDE}_{MILLIQUAS}_{REDSHIFT}_result.npy")

ts_mil_tur_xray = np.load(data_results_path / f"{TS_FILENAME}_{MILLIPEDE}_{TURIN}_{XRAY}.npy")
ts_mil_tur_xray_result = np.load(data_results_path / f"{TS_FILENAME}_{MILLIPEDE}_{TURIN}_{XRAY}_result.npy")

print(f"Plotting {SPLINEMPE} with {REDSHIFT}...")

plt.subplots(figsize=FIGSIZE)
logbins = np.logspace(0, np.log10(3e3), NBINS)
plt.hist(
    -ts_spl_tur_z,
    histtype=HISTTYPE,
    bins=logbins,
    linewidth = LINEWIDTH,
    label= TURIN_LABEL
)
plt.hist(
    -ts_spl_mil_z,
    histtype=HISTTYPE,
    bins=logbins,
    linewidth = LINEWIDTH,
    label= MILLIQUAS_LABEL
)
plt.axvline(
    -ts_spl_tur_z_result,
    color=AXVCOLOR,
    linestyle=AXVLINESTYLE,
    label=AXVLABEL
)
plt.gca().invert_xaxis()
plt.xlabel(XLABEL, fontsize=FONTSIZE)
plt.ylabel(YLABEL, fontsize=FONTSIZE)
plt.title(TITLE, fontsize=FONTSIZE)
plt.xscale(AXISSCALE)
plt.yscale(AXISSCALE)
plt.legend()

print("Saving plot...")

figname = f"{FIGNAME}_{SPLINEMPE}_{REDSHIFT}"
plt.savefig(figures_path / figname, bbox_inches = BBOX_INCHES)
figname_pdf = f"{figname}.pdf"
plt.savefig(figures_path / figname_pdf, bbox_inches = BBOX_INCHES)
plt.close()

print(f"Plotting {SPLINEMPE} with {XRAY}...")

plt.subplots(figsize=FIGSIZE)
logbins = np.logspace(0, np.log10(2e3), NBINS)
plt.hist(
    -ts_spl_tur_xray-np.min(-ts_spl_tur_xray)+1,
    histtype=HISTTYPE,
    bins=logbins,
    linewidth = LINEWIDTH,
    label= TURIN_LABEL
)
plt.axvline(
    -ts_spl_tur_xray_result-np.min(-ts_spl_tur_xray)+1,
    color=AXVCOLOR,
    linestyle=AXVLINESTYLE,
    label=AXVLABEL
)
plt.gca().invert_xaxis()
plt.xlabel(XLABEL, fontsize=FONTSIZE)
plt.ylabel(YLABEL, fontsize=FONTSIZE)
plt.title(XRAYTITLE, fontsize=FONTSIZE)
plt.xscale(AXISSCALE)
plt.yscale(AXISSCALE)
plt.legend()

print("Saving plot...")

figname = f"{FIGNAME}_{SPLINEMPE}_{XRAY}"
plt.savefig(figures_path / figname, bbox_inches = BBOX_INCHES)
figname_pdf = f"{figname}.pdf"
plt.savefig(figures_path / figname_pdf, bbox_inches = BBOX_INCHES)
plt.close()

print(f"Plotting {MILLIPEDE} with {REDSHIFT}...")

plt.subplots(figsize=FIGSIZE)
logbins = np.logspace(0, np.log10(60), NBINS)
plt.hist(
    -ts_mil_tur_z,
    histtype=HISTTYPE,
    bins=logbins,
    linewidth = LINEWIDTH,
    label= TURIN_LABEL
)
plt.hist(
    -ts_mil_mil_z,
    histtype=HISTTYPE,
    bins=logbins,
    linewidth = LINEWIDTH,
    label= MILLIQUAS_LABEL
)
plt.axvline(
    -ts_mil_tur_z_result,
    color=AXVCOLOR,
    linestyle=AXVLINESTYLE,
    label=AXVLABEL
)
plt.gca().invert_xaxis()
plt.xlabel(XLABEL, fontsize=FONTSIZE)
plt.ylabel(YLABEL, fontsize=FONTSIZE)
plt.title(MILLIPEDETITLE, fontsize=FONTSIZE)
plt.xscale(AXISSCALE)
plt.yscale(AXISSCALE)
plt.legend()

print("Saving plot...")

figname = f"{FIGNAME}_{MILLIPEDE}_{REDSHIFT}"
plt.savefig(figures_path / figname, bbox_inches = BBOX_INCHES)
figname_pdf = f"{figname}.pdf"
plt.savefig(figures_path / figname_pdf, bbox_inches = BBOX_INCHES)
plt.close()

print(f"Plotting {MILLIPEDE} with {XRAY}...")

plt.subplots(figsize=FIGSIZE)
logbins = np.logspace(0, np.log10(2e3), NBINS)
plt.hist(
    -ts_mil_tur_xray-np.min(-ts_mil_tur_xray)+1,
    histtype=HISTTYPE,
    bins=logbins,
    linewidth = LINEWIDTH,
    label= TURIN_LABEL
)
plt.axvline(
    -ts_mil_tur_xray_result-np.min(-ts_mil_tur_xray)+1,
    color=AXVCOLOR,
    linestyle=AXVLINESTYLE,
    label=AXVLABEL
)
plt.gca().invert_xaxis()
plt.xlabel(XLABEL, fontsize=FONTSIZE)
plt.ylabel(YLABEL, fontsize=FONTSIZE)
plt.title(MILLIPEDEXRAYTITLE, fontsize=FONTSIZE)
plt.xscale(AXISSCALE)
plt.yscale(AXISSCALE)
plt.legend()

print("Saving plot...")

figname = f"{FIGNAME}_{MILLIPEDE}_{XRAY}"
plt.savefig(figures_path / figname, bbox_inches = BBOX_INCHES)
figname_pdf = f"{figname}.pdf"
plt.savefig(figures_path / figname_pdf, bbox_inches = BBOX_INCHES)
plt.close()