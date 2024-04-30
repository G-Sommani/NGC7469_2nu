import numpy as np
import matplotlib.pyplot as plt
import config as cfg
from loading_functions import Loader

loader = Loader()
data_results_path = loader.data_results_path
figures_path = loader.figures_path

print("Loading results...")

ts_spl_tur_z = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.TURIN}_{cfg.REDSHIFT}.npy"
)
ts_spl_tur_z_result = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.TURIN}_{cfg.REDSHIFT}_result.npy"
)

ts_spl_tur_z_alt = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.TURIN}_{cfg.REDSHIFT}_alternative_hypothesis.npy"
)
ts_spl_tur_z_alt_result = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.TURIN}_{cfg.REDSHIFT}_alternative_hypothesis_result.npy"
)

ts_spl_mil_z = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.MILLIQUAS}_{cfg.REDSHIFT}.npy"
)
ts_spl_mil_z_result = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.MILLIQUAS}_{cfg.REDSHIFT}_result.npy"
)

ts_spl_tur_xray = np.load(
    data_results_path / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.TURIN}_{cfg.XRAY}.npy"
)
ts_spl_tur_xray_result = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.TURIN}_{cfg.XRAY}_result.npy"
)

ts_mil_tur_z = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.MILLIPEDE}_{cfg.TURIN}_{cfg.REDSHIFT}.npy"
)
ts_mil_tur_z_result = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.MILLIPEDE}_{cfg.TURIN}_{cfg.REDSHIFT}_result.npy"
)

ts_mil_mil_z = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.MILLIPEDE}_{cfg.MILLIQUAS}_{cfg.REDSHIFT}.npy"
)
ts_mil_mil_z_result = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.MILLIPEDE}_{cfg.MILLIQUAS}_{cfg.REDSHIFT}_result.npy"
)

ts_mil_tur_xray = np.load(
    data_results_path / f"{cfg.TS_FILENAME}_{cfg.MILLIPEDE}_{cfg.TURIN}_{cfg.XRAY}.npy"
)
ts_mil_tur_xray_result = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.MILLIPEDE}_{cfg.TURIN}_{cfg.XRAY}_result.npy"
)

print(f"Plotting {cfg.SPLINEMPE} with {cfg.REDSHIFT}...")

plt.subplots(figsize=cfg.FIGSIZE_TS)
logbins = list(np.logspace(0, np.log10(3e3), cfg.NBINS_TS))
plt.hist(
    -ts_spl_tur_z,
    histtype=cfg.HISTTYPE_TS,
    bins=logbins,
    linewidth=cfg.LINEWIDTH_TS,
    label=cfg.TURIN_LABEL_TS,
)
plt.hist(
    -ts_spl_mil_z,
    histtype=cfg.HISTTYPE_TS,
    bins=logbins,
    linewidth=cfg.LINEWIDTH_TS,
    label=cfg.MILLIQUAS_LABEL_TS,
)
plt.axvline(
    -ts_spl_tur_z_result,
    color=cfg.AXVCOLOR_TS,
    linestyle=cfg.AXVLINESTYLE_TS,
    label=cfg.AXVLABEL_TS,
)
plt.gca().invert_xaxis()
plt.xlabel(cfg.XLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.ylabel(cfg.YLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.title(cfg.TITLE_TS, fontsize=cfg.FONTSIZE_TS)
plt.xscale(cfg.AXISSCALE_TS)
plt.yscale(cfg.AXISSCALE_TS)
plt.legend()

print("Saving plot...")

figname = f"{cfg.FIGNAME_TS}_{cfg.SPLINEMPE}_{cfg.REDSHIFT}"
plt.savefig(figures_path / figname, bbox_inches=cfg.BBOX_INCHES)
figname_pdf = f"{figname}.pdf"
plt.savefig(figures_path / figname_pdf, bbox_inches=cfg.BBOX_INCHES)
plt.close()

print(f"Plotting {cfg.SPLINEMPE} with {cfg.XRAY}...")

plt.subplots(figsize=cfg.FIGSIZE_TS)
logbins = list(np.logspace(0, np.log10(2e3), cfg.NBINS_TS))
plt.hist(
    -ts_spl_tur_xray - np.min(-ts_spl_tur_xray) + 1,
    histtype=cfg.HISTTYPE_TS,
    bins=logbins,
    linewidth=cfg.LINEWIDTH_TS,
    label=cfg.TURIN_LABEL_TS,
)
plt.axvline(
    -ts_spl_tur_xray_result - np.min(-ts_spl_tur_xray) + 1,
    color=cfg.AXVCOLOR_TS,
    linestyle=cfg.AXVLINESTYLE_TS,
    label=cfg.AXVLABEL_TS,
)
plt.gca().invert_xaxis()
plt.xlabel(cfg.XLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.ylabel(cfg.YLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.title(cfg.XRAYTITLE_TS, fontsize=cfg.FONTSIZE_TS)
plt.xscale(cfg.AXISSCALE_TS)
plt.yscale(cfg.AXISSCALE_TS)
plt.legend()

print("Saving plot...")

figname = f"{cfg.FIGNAME_TS}_{cfg.SPLINEMPE}_{cfg.XRAY}"
plt.savefig(figures_path / figname, bbox_inches=cfg.BBOX_INCHES)
figname_pdf = f"{figname}.pdf"
plt.savefig(figures_path / figname_pdf, bbox_inches=cfg.BBOX_INCHES)
plt.close()

print(f"Plotting {cfg.MILLIPEDE} with {cfg.REDSHIFT}...")

plt.subplots(figsize=cfg.FIGSIZE_TS)
logbins = list(np.logspace(0, np.log10(60), cfg.NBINS_TS))
plt.hist(
    -ts_mil_tur_z,
    histtype=cfg.HISTTYPE_TS,
    bins=logbins,
    linewidth=cfg.LINEWIDTH_TS,
    label=cfg.TURIN_LABEL_TS,
)
plt.hist(
    -ts_mil_mil_z,
    histtype=cfg.HISTTYPE_TS,
    bins=logbins,
    linewidth=cfg.LINEWIDTH_TS,
    label=cfg.MILLIQUAS_LABEL_TS,
)
plt.axvline(
    -ts_mil_tur_z_result,
    color=cfg.AXVCOLOR_TS,
    linestyle=cfg.AXVLINESTYLE_TS,
    label=cfg.AXVLABEL_TS,
)
plt.gca().invert_xaxis()
plt.xlabel(cfg.XLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.ylabel(cfg.YLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.title(cfg.MILLIPEDETITLE_TS, fontsize=cfg.FONTSIZE_TS)
plt.xscale(cfg.AXISSCALE_TS)
plt.yscale(cfg.AXISSCALE_TS)
plt.legend()

print("Saving plot...")

figname = f"{cfg.FIGNAME_TS}_{cfg.MILLIPEDE}_{cfg.REDSHIFT}"
plt.savefig(figures_path / figname, bbox_inches=cfg.BBOX_INCHES)
figname_pdf = f"{figname}.pdf"
plt.savefig(figures_path / figname_pdf, bbox_inches=cfg.BBOX_INCHES)
plt.close()

print(f"Plotting {cfg.MILLIPEDE} with {cfg.XRAY}...")

plt.subplots(figsize=cfg.FIGSIZE_TS)
logbins = list(np.logspace(0, np.log10(2e3), cfg.NBINS_TS))
plt.hist(
    -ts_mil_tur_xray - np.min(-ts_mil_tur_xray) + 1,
    histtype=cfg.HISTTYPE_TS,
    bins=logbins,
    linewidth=cfg.LINEWIDTH_TS,
    label=cfg.TURIN_LABEL_TS,
)
plt.axvline(
    -ts_mil_tur_xray_result - np.min(-ts_mil_tur_xray) + 1,
    color=cfg.AXVCOLOR_TS,
    linestyle=cfg.AXVLINESTYLE_TS,
    label=cfg.AXVLABEL_TS,
)
plt.gca().invert_xaxis()
plt.xlabel(cfg.XLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.ylabel(cfg.YLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.title(cfg.MILLIPEDEXRAYTITLE_TS, fontsize=cfg.FONTSIZE_TS)
plt.xscale(cfg.AXISSCALE_TS)
plt.yscale(cfg.AXISSCALE_TS)
plt.legend()

print("Saving plot...")

figname = f"{cfg.FIGNAME_TS}_{cfg.MILLIPEDE}_{cfg.XRAY}"
plt.savefig(figures_path / figname, bbox_inches=cfg.BBOX_INCHES)
figname_pdf = f"{figname}.pdf"
plt.savefig(figures_path / figname_pdf, bbox_inches=cfg.BBOX_INCHES)
plt.close()

print(
    f"Plotting alternative hypothesis vs background hypothesis for {cfg.SPLINEMPE} with {cfg.TURIN} and {cfg.REDSHIFT}..."
)

plt.subplots(figsize=cfg.FIGSIZE_TS)
logbins = list(np.logspace(0, np.log10(3e3), cfg.NBINS_TS))
plt.hist(
    -ts_spl_tur_z,
    histtype=cfg.HISTTYPE_TS,
    bins=logbins,
    linewidth=cfg.LINEWIDTH_TS,
    label=f"{cfg.TURIN_LABEL_TS} background",
)
plt.hist(
    -ts_spl_tur_z_alt,
    histtype=cfg.HISTTYPE_TS,
    bins=logbins,
    linewidth=cfg.LINEWIDTH_TS,
    label=f"{cfg.TURIN_LABEL_TS} alternative hypothesis",
)
plt.axvline(
    -ts_spl_tur_z_result,
    color=cfg.AXVCOLOR_TS,
    linestyle=cfg.AXVLINESTYLE_TS,
    label=cfg.AXVLABEL_TS,
)
plt.gca().invert_xaxis()
plt.xlabel(cfg.XLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.ylabel(cfg.YLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.title(cfg.TITLE_TS, fontsize=cfg.FONTSIZE_TS)
plt.xscale(cfg.AXISSCALE_TS)
plt.yscale(cfg.AXISSCALE_TS)
plt.legend()

print("Saving plot...")

figname = f"{cfg.FIGNAME_TS}_{cfg.SPLINEMPE}_{cfg.TURIN}_{cfg.REDSHIFT}_background_vs_alternative_hypothesis"
plt.savefig(figures_path / figname, bbox_inches=cfg.BBOX_INCHES)
figname_pdf = f"{figname}.pdf"
plt.savefig(figures_path / figname_pdf, bbox_inches=cfg.BBOX_INCHES)
plt.close()
