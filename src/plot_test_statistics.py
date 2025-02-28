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
    / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.TURIN}_{cfg.REDSHIFT}_Background_hypothesis.npy"
)
ts_spl_tur_z_result = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.TURIN}_{cfg.REDSHIFT}_Background_hypothesis_result.npy"
)

ts_spl_tur_z_jitter = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.TURIN}_{cfg.REDSHIFT}_dec_jitter_Background_hypothesis.npy"
)
ts_spl_tur_z_jitter_result = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.TURIN}_{cfg.REDSHIFT}_dec_jitter_Background_hypothesis_result.npy"
)


ts_spl_tur_z_alt = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.TURIN}_{cfg.REDSHIFT}_Doublet_hypothesis.npy"
)
ts_spl_tur_z_alt_result = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.TURIN}_{cfg.REDSHIFT}_Doublet_hypothesis_result.npy"
)

ts_spl_tur_z_alt_sing_inj = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.TURIN}_{cfg.REDSHIFT}_Singlet-inj_hypothesis.npy"
)
ts_spl_tur_z_alt_sing_inj_result = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.TURIN}_{cfg.REDSHIFT}_Singlet-inj_hypothesis_result.npy"
)

ts_spl_tur_z_alt_inj = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.TURIN}_{cfg.REDSHIFT}_Doublet-inj_hypothesis.npy"
)
ts_spl_tur_z_alt_inj_result = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.TURIN}_{cfg.REDSHIFT}_Doublet-inj_hypothesis_result.npy"
)
zs_tur = np.load(
    data_results_path / f"{cfg.REDSHIFT}_sources_{cfg.TURIN}_Doublet_hypothesis.npy"
)

ts_spl_mil_z = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.MILLIQUAS}_{cfg.REDSHIFT}_Background_hypothesis.npy"
)
ts_spl_mil_z_result = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.MILLIQUAS}_{cfg.REDSHIFT}_Background_hypothesis_result.npy"
)

ts_spl_mil_z_alt = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.MILLIQUAS}_{cfg.REDSHIFT}_Doublet_hypothesis.npy"
)
ts_spl_mil_z_alt_result = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.MILLIQUAS}_{cfg.REDSHIFT}_Doublet_hypothesis_result.npy"
)

ts_spl_tur_xray = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.TURIN}_{cfg.XRAY}_Background_hypothesis.npy"
)
ts_spl_tur_xray_result = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.TURIN}_{cfg.XRAY}_Background_hypothesis_result.npy"
)

ts_spl_tur_xray_alt = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.TURIN}_{cfg.XRAY}_Doublet_hypothesis.npy"
)
ts_spl_tur_xray_alt_result = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.SPLINEMPE}_{cfg.TURIN}_{cfg.XRAY}_Doublet_hypothesis_result.npy"
)

ts_mil_tur_z = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.MILLIPEDE}_{cfg.TURIN}_{cfg.REDSHIFT}_Background_hypothesis.npy"
)
ts_mil_tur_z_result = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.MILLIPEDE}_{cfg.TURIN}_{cfg.REDSHIFT}_Background_hypothesis_result.npy"
)

ts_mil_tur_z_alt = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.MILLIPEDE}_{cfg.TURIN}_{cfg.REDSHIFT}_Doublet_hypothesis.npy"
)
ts_mil_tur_z_alt_result = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.MILLIPEDE}_{cfg.TURIN}_{cfg.REDSHIFT}_Doublet_hypothesis_result.npy"
)

ts_mil_mil_z = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.MILLIPEDE}_{cfg.MILLIQUAS}_{cfg.REDSHIFT}_Background_hypothesis.npy"
)
ts_mil_mil_z_result = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.MILLIPEDE}_{cfg.MILLIQUAS}_{cfg.REDSHIFT}_Background_hypothesis_result.npy"
)

ts_mil_mil_z_alt = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.MILLIPEDE}_{cfg.MILLIQUAS}_{cfg.REDSHIFT}_Doublet_hypothesis.npy"
)
ts_mil_mil_z_alt_result = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.MILLIPEDE}_{cfg.MILLIQUAS}_{cfg.REDSHIFT}_Doublet_hypothesis_result.npy"
)

ts_mil_tur_xray = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.MILLIPEDE}_{cfg.TURIN}_{cfg.XRAY}_Background_hypothesis.npy"
)
ts_mil_tur_xray_result = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.MILLIPEDE}_{cfg.TURIN}_{cfg.XRAY}_Background_hypothesis_result.npy"
)

ts_mil_tur_xray_alt = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.MILLIPEDE}_{cfg.TURIN}_{cfg.XRAY}_Doublet_hypothesis.npy"
)
ts_mil_tur_xray_alt_result = np.load(
    data_results_path
    / f"{cfg.TS_FILENAME}_{cfg.MILLIPEDE}_{cfg.TURIN}_{cfg.XRAY}_Doublet_hypothesis_result.npy"
)

print("\n\n*** Background Hypotheses ***\n\n")
print(f"Plotting {cfg.SPLINEMPE} with {cfg.REDSHIFT}...")

plt.subplots(figsize=cfg.FIGSIZE_TS)
logbins = list(np.logspace(0, np.log10(7e2), cfg.NBINS_TS))
plt.hist(
    -ts_spl_tur_z - np.min(-ts_spl_tur_z) + 1,
    histtype=cfg.HISTTYPE_TS,
    bins=logbins,
    linewidth=cfg.LINEWIDTH_TS,
    label=cfg.TURIN_LABEL_TS,
)
plt.hist(
    -ts_spl_mil_z - np.min(-ts_spl_tur_z) + 1,
    histtype=cfg.HISTTYPE_TS,
    bins=logbins,
    linewidth=cfg.LINEWIDTH_TS,
    label=cfg.MILLIQUAS_LABEL_TS,
)
plt.axvline(
    -ts_spl_tur_z_result - np.min(-ts_spl_tur_z) + 1,
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
plt.savefig(figures_path / figname_pdf, bbox_inches=cfg.BBOX_INCHES, dpi=200)
plt.close()

print(f"Plotting {cfg.SPLINEMPE} with {cfg.REDSHIFT} (declination jitter)...")

plt.subplots(figsize=cfg.FIGSIZE_TS)
logbins = list(np.logspace(0, np.log10(4e2), cfg.NBINS_TS))
plt.hist(
    -ts_spl_tur_z_jitter,
    histtype=cfg.HISTTYPE_TS,
    bins=logbins,
    linewidth=cfg.LINEWIDTH_TS,
    label=cfg.TURIN_LABEL_TS,
)
plt.hist(
    -ts_spl_tur_z,
    histtype=cfg.HISTTYPE_TS,
    bins=logbins,
    linewidth=cfg.LINEWIDTH_TS,
    label=cfg.TURIN_LABEL_TS,
)
plt.axvline(
    -ts_spl_tur_z_jitter_result,
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

figname = f"{cfg.FIGNAME_TS}_{cfg.SPLINEMPE}_{cfg.REDSHIFT}_dec_jitter"
plt.savefig(figures_path / figname, bbox_inches=cfg.BBOX_INCHES)
figname_pdf = f"{figname}.pdf"
plt.savefig(figures_path / figname_pdf, bbox_inches=cfg.BBOX_INCHES, dpi=200)
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

print("\n\n*** Doublet Hypotheses ***\n\n")

print(f"Plotting {cfg.SPLINEMPE} with {cfg.TURIN} and {cfg.REDSHIFT}...")

plt.subplots(figsize=cfg.FIGSIZE_TS)
logbins = list(np.logspace(0, np.log10(4e2), cfg.NBINS_TS))
alpha = 0.7
plt.hist(
    -ts_spl_tur_z_alt - np.min(-ts_spl_tur_z_alt) + 1,
    histtype=cfg.HISTTYPE_TS_ALT,
    bins=logbins,
    alpha=alpha,
    edgecolor="blue",
    linewidth=cfg.LINEWIDTH_TS,
    label=f"Alternative hypothesis distribution\nTurin catalog",
)
plt.hist(
    -ts_spl_tur_z - np.min(-ts_spl_tur_z_alt) + 1,
    # - (ts_spl_tur_z_alt_result - ts_spl_tur_z_result)
    histtype=cfg.HISTTYPE_TS_ALT,
    bins=logbins,
    alpha=alpha,
    edgecolor="red",
    linewidth=cfg.LINEWIDTH_TS,
    label=f"{cfg.TURIN_LABEL_TS}",
)
plt.axvline(
    -ts_spl_tur_z_alt_result - np.min(-ts_spl_tur_z_alt) + 1,
    color=cfg.AXVCOLOR_TS,
    linestyle=cfg.AXVLINESTYLE_TS,
    label=cfg.AXVLABEL_TS,
)
plt.gca().invert_xaxis()
plt.xlabel(cfg.XLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.ylabel(cfg.YLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.title(f"{cfg.TITLE_UNBIASEDNESS}", fontsize=cfg.FONTSIZE_TS)
plt.xscale(cfg.AXISSCALE_TS)
plt.yscale(cfg.AXISSCALE_TS)
plt.legend()

print("Saving plot...")

figname = f"{cfg.FIGNAME_TS}_{cfg.SPLINEMPE}_{cfg.TURIN}_{cfg.REDSHIFT}_background_vs_alternative_hypothesis_doublet"
plt.savefig(figures_path / figname, bbox_inches=cfg.BBOX_INCHES)
figname_pdf = f"{figname}.pdf"
plt.savefig(figures_path / figname_pdf, bbox_inches=cfg.BBOX_INCHES)
plt.close()

print(
    f"Plotting {cfg.SPLINEMPE} with {cfg.TURIN} and {cfg.REDSHIFT} (injected case)..."
)

plt.subplots(figsize=cfg.FIGSIZE_TS)
logbins = list(np.logspace(0, np.log10(4e2), cfg.NBINS_TS))
alpha = 0.7
plt.hist(
    -ts_spl_tur_z_alt_inj - np.min(-ts_spl_tur_z_alt_inj) + 1,
    histtype=cfg.HISTTYPE_TS_ALT,
    bins=logbins,
    alpha=alpha,
    edgecolor="blue",
    linewidth=cfg.LINEWIDTH_TS,
    label=f"Doublet injection",
)
plt.hist(
    -ts_spl_tur_z
    - np.min(-ts_spl_tur_z_alt_inj)
    + 1
    - (ts_spl_tur_z_alt_inj_result - ts_spl_tur_z_result),
    histtype=cfg.HISTTYPE_TS_ALT,
    bins=logbins,
    alpha=alpha,
    edgecolor="red",
    linewidth=cfg.LINEWIDTH_TS,
    label=f"Null hypothesis distribution",
)
plt.axvline(
    -ts_spl_tur_z_alt_inj_result - np.min(-ts_spl_tur_z_alt_inj) + 1,
    color=cfg.AXVCOLOR_TS,
    linestyle=cfg.AXVLINESTYLE_TS,
    label=cfg.AXVLABEL_TS,
)
plt.gca().invert_xaxis()
plt.xlabel(cfg.XLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.ylabel(cfg.YLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.title(
    f"Sensitivity to the injection of doublet coincidences", fontsize=cfg.FONTSIZE_TS
)
plt.xscale(cfg.AXISSCALE_TS)
plt.yscale(cfg.AXISSCALE_TS)
plt.legend()

print("Saving plot...")

figname = f"{cfg.FIGNAME_TS}_{cfg.SPLINEMPE}_{cfg.TURIN}_{cfg.REDSHIFT}_background_vs_alternative_hypothesis_doublet_injected"
plt.savefig(figures_path / figname, bbox_inches=cfg.BBOX_INCHES)
figname_pdf = f"{figname}.pdf"
plt.savefig(figures_path / figname_pdf, bbox_inches=cfg.BBOX_INCHES)
plt.close()

print(
    f"Plotting {cfg.SPLINEMPE} with {cfg.TURIN} and {cfg.REDSHIFT} (injected case with singlet)..."
)

plt.subplots(figsize=cfg.FIGSIZE_TS)
logbins = list(np.logspace(0, np.log10(4e2), cfg.NBINS_TS))
alpha = 1
plt.hist(
    -ts_spl_tur_z
    - np.min(-ts_spl_tur_z_alt_inj)
    + 1
    - (ts_spl_tur_z_alt_inj_result - ts_spl_tur_z_result),
    histtype=cfg.HISTTYPE_TS_ALT,
    bins=logbins,
    alpha=0.5,
    edgecolor="black",
    color="grey",
    linewidth=cfg.LINEWIDTH_TS - 1,
    label=f"Null hypothesis",
)
plt.hist(
    -ts_spl_tur_z_alt_sing_inj - np.min(-ts_spl_tur_z_alt_inj) + 1,
    histtype=cfg.HISTTYPE_TS,
    bins=logbins,
    alpha=alpha,
    edgecolor="tab:blue",
    color="tab:green",
    linewidth=cfg.LINEWIDTH_TS,
    label=f"Singlet injection",
)
plt.hist(
    -ts_spl_tur_z_alt_inj - np.min(-ts_spl_tur_z_alt_inj) + 1,
    histtype=cfg.HISTTYPE_TS,
    bins=logbins,
    alpha=alpha,
    edgecolor="tab:orange",
    linewidth=cfg.LINEWIDTH_TS,
    label=f"Doublet injection",
)
plt.axvline(
    -ts_spl_tur_z_alt_sing_inj_result - np.min(-ts_spl_tur_z_alt_inj) + 1,
    color="red",
    linestyle="dotted",
    linewidth=2.5,
    label=cfg.AXVLABEL_TS,
)
plt.gca().invert_xaxis()
plt.xlabel(cfg.XLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.ylabel(cfg.YLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.title(
    f"Sensitivity to the injection\nof doublet and singlet coincidences",
    fontsize=cfg.FONTSIZE_TS,
)
plt.xscale(cfg.AXISSCALE_TS)
plt.yscale(cfg.AXISSCALE_TS)
plt.legend()

print("Saving plot...")

figname = f"{cfg.FIGNAME_TS}_{cfg.SPLINEMPE}_{cfg.TURIN}_{cfg.REDSHIFT}_background_vs_alternative_hypothesis_singlet_injected"
plt.savefig(figures_path / figname, bbox_inches=cfg.BBOX_INCHES)
figname_pdf = f"{figname}.pdf"
plt.savefig(figures_path / figname_pdf, bbox_inches=cfg.BBOX_INCHES)
plt.close()

print(f"Plotting {cfg.REDSHIFT} of selected sources under {cfg.TURIN} catalog...")

plt.figure(figsize=cfg.FIGSIZE_TS)
logbins = list(np.logspace(np.log10(min(zs_tur)), np.log10(max(zs_tur)), 17))
plt.hist(
    zs_tur,
    bins=logbins,
    linewidth=2,
    edgecolor="blue",
    alpha=0.7,
    histtype="stepfilled",
    label="Selected redshifts",
)
plt.axvline(0.016, color="red", linestyle="--", linewidth=2, label="NGC 7469")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Redshift", fontsize="large")
plt.ylabel("Number of selections", fontsize="large")
plt.title("Selected resdhifts", fontsize="large")
plt.legend()

print("Saving plot...")

plt.savefig(figures_path / "selected_redshifts_turin", bbox_inches=cfg.BBOX_INCHES)
plt.savefig(figures_path / "selected_redshifts_turin.pdf", bbox_inches=cfg.BBOX_INCHES)

print(f"Plotting {cfg.SPLINEMPE} with {cfg.TURIN} and {cfg.XRAY}...")

plt.subplots(figsize=cfg.FIGSIZE_TS)
logbins = list(np.logspace(0, np.log10(2e3), cfg.NBINS_TS))
alpha = 0.7
plt.hist(
    -ts_spl_tur_xray_alt - np.min(-ts_spl_tur_xray_alt) + 1,
    histtype=cfg.HISTTYPE_TS_ALT,
    bins=logbins,
    alpha=alpha,
    edgecolor="blue",
    linewidth=cfg.LINEWIDTH_TS,
    label=f"Alternative hypothesis distribution\nTurin catalog",
)
plt.hist(
    -ts_spl_tur_xray
    - np.min(-ts_spl_tur_xray_alt)
    - (ts_spl_tur_xray_alt_result - ts_spl_tur_xray_result)
    + 1,
    histtype=cfg.HISTTYPE_TS_ALT,
    bins=logbins,
    alpha=alpha,
    edgecolor="red",
    linewidth=cfg.LINEWIDTH_TS,
    label=f"{cfg.TURIN_LABEL_TS}",
)
plt.axvline(
    -ts_spl_tur_xray_alt_result - np.min(-ts_spl_tur_xray_alt) + 1,
    color=cfg.AXVCOLOR_TS,
    linestyle=cfg.AXVLINESTYLE_TS,
    label=cfg.AXVLABEL_TS,
)
plt.gca().invert_xaxis()
plt.xlabel(cfg.XLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.ylabel(cfg.YLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.title(f"{cfg.TITLE_UNBIASEDNESS}", fontsize=cfg.FONTSIZE_TS)
plt.xscale(cfg.AXISSCALE_TS)
plt.yscale(cfg.AXISSCALE_TS)
plt.legend()

print("Saving plot...")

figname = f"{cfg.FIGNAME_TS}_{cfg.SPLINEMPE}_{cfg.TURIN}_{cfg.XRAY}_background_vs_alternative_hypothesis_doublet"
plt.savefig(figures_path / figname, bbox_inches=cfg.BBOX_INCHES)
figname_pdf = f"{figname}.pdf"
plt.savefig(figures_path / figname_pdf, bbox_inches=cfg.BBOX_INCHES)
plt.close()

print(f"Plotting {cfg.SPLINEMPE} with {cfg.MILLIQUAS} and {cfg.REDSHIFT}...")

plt.subplots(figsize=cfg.FIGSIZE_TS)
logbins = list(np.logspace(0, np.log10(2e3), cfg.NBINS_TS))
alpha = 0.7
plt.hist(
    -ts_spl_mil_z_alt - np.min(-ts_spl_mil_z_alt) + 1,
    histtype=cfg.HISTTYPE_TS_ALT,
    bins=logbins,
    alpha=alpha,
    edgecolor="blue",
    linewidth=cfg.LINEWIDTH_TS,
    label=f"Alternative hypothesis distribution\nMilliquas catalog",
)
plt.hist(
    -ts_spl_mil_z
    - np.min(-ts_spl_mil_z_alt)
    - (ts_spl_mil_z_alt_result - ts_spl_mil_z_result)
    + 1,
    histtype=cfg.HISTTYPE_TS_ALT,
    bins=logbins,
    alpha=alpha,
    edgecolor="red",
    linewidth=cfg.LINEWIDTH_TS,
    label=f"{cfg.MILLIQUAS_LABEL_TS}",
)
plt.axvline(
    -ts_spl_mil_z_alt_result - np.min(-ts_spl_mil_z_alt) + 1,
    color=cfg.AXVCOLOR_TS,
    linestyle=cfg.AXVLINESTYLE_TS,
    label=cfg.AXVLABEL_TS,
)
plt.gca().invert_xaxis()
plt.xlabel(cfg.XLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.ylabel(cfg.YLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.title(f"{cfg.TITLE_UNBIASEDNESS}", fontsize=cfg.FONTSIZE_TS)
plt.xscale(cfg.AXISSCALE_TS)
plt.yscale(cfg.AXISSCALE_TS)
plt.legend()

print("Saving plot...")

figname = f"{cfg.FIGNAME_TS}_{cfg.SPLINEMPE}_{cfg.MILLIQUAS}_{cfg.REDSHIFT}_background_vs_alternative_hypothesis_doublet"
plt.savefig(figures_path / figname, bbox_inches=cfg.BBOX_INCHES)
figname_pdf = f"{figname}.pdf"
plt.savefig(figures_path / figname_pdf, bbox_inches=cfg.BBOX_INCHES)
plt.close()


print(f"Plotting {cfg.MILLIPEDE} with {cfg.TURIN} and {cfg.REDSHIFT}...")

plt.subplots(figsize=cfg.FIGSIZE_TS)
logbins = list(np.logspace(0, np.log10(2e3), cfg.NBINS_TS))
alpha = 0.7
plt.hist(
    -ts_mil_tur_z_alt - np.min(-ts_mil_tur_z_alt) + 1,
    histtype=cfg.HISTTYPE_TS_ALT,
    bins=logbins,
    alpha=alpha,
    edgecolor="blue",
    linewidth=cfg.LINEWIDTH_TS,
    label=f"Alternative hypothesis distribution\nTurin catalog",
)
plt.hist(
    -ts_mil_tur_z
    - np.min(-ts_mil_tur_z_alt)
    - (ts_mil_tur_z_alt_result - ts_mil_tur_z_result)
    + 1,
    histtype=cfg.HISTTYPE_TS_ALT,
    bins=logbins,
    alpha=alpha,
    edgecolor="red",
    linewidth=cfg.LINEWIDTH_TS,
    label=f"{cfg.TURIN_LABEL_TS}",
)
plt.axvline(
    -ts_mil_tur_z_alt_result - np.min(-ts_mil_tur_z_alt) + 1,
    color=cfg.AXVCOLOR_TS,
    linestyle=cfg.AXVLINESTYLE_TS,
    label=cfg.AXVLABEL_TS,
)
plt.gca().invert_xaxis()
plt.xlabel(cfg.XLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.ylabel(cfg.YLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.title(f"{cfg.TITLE_UNBIASEDNESS}", fontsize=cfg.FONTSIZE_TS)
plt.xscale(cfg.AXISSCALE_TS)
plt.yscale(cfg.AXISSCALE_TS)
plt.legend()

print("Saving plot...")

figname = f"{cfg.FIGNAME_TS}_{cfg.MILLIPEDE}_{cfg.TURIN}_{cfg.REDSHIFT}_background_vs_alternative_hypothesis_doublet"
plt.savefig(figures_path / figname, bbox_inches=cfg.BBOX_INCHES)
figname_pdf = f"{figname}.pdf"
plt.savefig(figures_path / figname_pdf, bbox_inches=cfg.BBOX_INCHES)
plt.close()


print(f"Plotting {cfg.MILLIPEDE} with {cfg.TURIN} and {cfg.XRAY}...")

plt.subplots(figsize=cfg.FIGSIZE_TS)
logbins = list(np.logspace(0, np.log10(2e3), cfg.NBINS_TS))
alpha = 0.7
plt.hist(
    -ts_mil_tur_xray_alt - np.nanmin(-ts_mil_tur_xray_alt) + 1,
    histtype=cfg.HISTTYPE_TS_ALT,
    bins=logbins,
    alpha=alpha,
    edgecolor="blue",
    linewidth=cfg.LINEWIDTH_TS,
    label=f"Alternative hypothesis distribution\nTurin catalog",
)
plt.hist(
    -ts_mil_tur_xray
    - np.nanmin(-ts_mil_tur_xray_alt)
    - (ts_mil_tur_xray_alt_result - ts_mil_tur_xray_result)
    + 1,
    histtype=cfg.HISTTYPE_TS_ALT,
    bins=logbins,
    alpha=alpha,
    edgecolor="red",
    linewidth=cfg.LINEWIDTH_TS,
    label=f"{cfg.TURIN_LABEL_TS}",
)
plt.axvline(
    -ts_mil_tur_xray_alt_result - np.nanmin(-ts_mil_tur_xray_alt) + 1,
    color=cfg.AXVCOLOR_TS,
    linestyle=cfg.AXVLINESTYLE_TS,
    label=cfg.AXVLABEL_TS,
)
plt.gca().invert_xaxis()
plt.xlabel(cfg.XLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.ylabel(cfg.YLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.title(f"{cfg.TITLE_UNBIASEDNESS}", fontsize=cfg.FONTSIZE_TS)
plt.xscale(cfg.AXISSCALE_TS)
plt.yscale(cfg.AXISSCALE_TS)
plt.legend()

print("Saving plot...")

figname = f"{cfg.FIGNAME_TS}_{cfg.MILLIPEDE}_{cfg.TURIN}_{cfg.XRAY}_background_vs_alternative_hypothesis_doublet"
plt.savefig(figures_path / figname, bbox_inches=cfg.BBOX_INCHES)
figname_pdf = f"{figname}.pdf"
plt.savefig(figures_path / figname_pdf, bbox_inches=cfg.BBOX_INCHES)
plt.close()


print(f"Plotting {cfg.MILLIPEDE} with {cfg.MILLIQUAS} and {cfg.REDSHIFT}...")

plt.subplots(figsize=cfg.FIGSIZE_TS)
logbins = list(np.logspace(0, np.log10(2e3), cfg.NBINS_TS))
alpha = 0.7
plt.hist(
    -ts_mil_mil_z_alt - np.min(-ts_mil_mil_z_alt) + 1,
    histtype=cfg.HISTTYPE_TS_ALT,
    bins=logbins,
    alpha=alpha,
    edgecolor="blue",
    linewidth=cfg.LINEWIDTH_TS,
    label=f"Alternative hypothesis distribution\nMilliquas catalog",
)
plt.hist(
    -ts_mil_mil_z
    - np.min(-ts_mil_mil_z_alt)
    - (ts_mil_mil_z_alt_result - ts_mil_mil_z_result)
    + 1,
    histtype=cfg.HISTTYPE_TS_ALT,
    bins=logbins,
    alpha=alpha,
    edgecolor="red",
    linewidth=cfg.LINEWIDTH_TS,
    label=f"{cfg.MILLIQUAS_LABEL_TS}",
)
plt.axvline(
    -ts_mil_mil_z_alt_result - np.min(-ts_mil_mil_z_alt) + 1,
    color=cfg.AXVCOLOR_TS,
    linestyle=cfg.AXVLINESTYLE_TS,
    label=cfg.AXVLABEL_TS,
)
plt.gca().invert_xaxis()
plt.xlabel(cfg.XLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.ylabel(cfg.YLABEL_TS, fontsize=cfg.FONTSIZE_TS)
plt.title(f"{cfg.TITLE_UNBIASEDNESS}", fontsize=cfg.FONTSIZE_TS)
plt.xscale(cfg.AXISSCALE_TS)
plt.yscale(cfg.AXISSCALE_TS)
plt.legend()

print("Saving plot...")

figname = f"{cfg.FIGNAME_TS}_{cfg.MILLIPEDE}_{cfg.MILLIQUAS}_{cfg.REDSHIFT}_background_vs_alternative_hypothesis_doublet"
plt.savefig(figures_path / figname, bbox_inches=cfg.BBOX_INCHES)
figname_pdf = f"{figname}.pdf"
plt.savefig(figures_path / figname_pdf, bbox_inches=cfg.BBOX_INCHES)
plt.close()
