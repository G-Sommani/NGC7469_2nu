import numpy as np
import os
from pathlib import Path
import matplotlib.pyplot as plt
import config as cfg

print("Definition of paths...")

# Definition of paths
cwd = Path(os.getcwd())
data_path = cwd / "../data"
figures_path = cwd / "../figures"

effective_area = np.genfromtxt(data_path / cfg.EFFECTIVE_AREA_FILENAME)
energy_bins = effective_area[:, cfg.EFFECTIVE_AREA_ENERGY_BINS_INDEX]
effective_area_array = np.array(
    [
        effective_area[:, cfg.EFFECTIVE_AREA_30_90_DEG_INDEX],
        effective_area[:, cfg.EFFECTIVE_AREA_0_30_DEG_INDEX],
        effective_area[:, cfg.EFFECTIVE_AREA_MIN5_0_DEG_INDEX],
        effective_area[:, cfg.EFFECTIVE_AREA_MIN30_MIN5_DEG_INDEX],
        effective_area[:, cfg.EFFECTIVE_AREA_MIN90_MIN30_DEG_INDEX],
    ]
)


def energy_factor(bin_index):
    """
    Estimate energy factor for expected number of detected neutrinos
    """
    e_max = energy_bins[bin_index + 1]
    e_min = energy_bins[bin_index]
    factor = (e_max - e_min) / (e_max * e_min)
    return factor


def area_energy_factor_calculator(a_index):
    """
    Estimate area and energy factor for expected number of detected neutrinos
    """
    eff_area = effective_area_array[a_index]
    factor = 0
    for k in range(len(energy_bins) - 1):
        element = eff_area[k] * energy_factor(k)
        factor += element
    return factor  # units: m^2 GeV^-1


area_energy_factors = np.array([])
for i in range(len(effective_area_array)):
    area_energy_factors = np.append(
        area_energy_factors, area_energy_factor_calculator(i)
    )
hubble_in_s = cfg.HUBBLE_CONSTANT / cfg.MPC_TO_KM
days = cfg.MJD_04102023 - cfg.MJD_GOLDBRONZE_START
seconds = (days / cfg.DAYS_IN_YEAR) * cfg.SECONDS_IN_YEAR


def expected_nu_from_source(z, dec, nu_flux=cfg.FLUX_NU):
    """
    Given the redshift and the declination of a source, determines the total
    number of expected neutrinos from the source
    """
    area_energy_factor = None
    if 90 >= dec > 30:
        area_energy_factor = area_energy_factors[cfg.EFFECTIVE_AREA_30_90_DEG_INDEX - 1]
    elif dec <= 30 and dec > 0:
        area_energy_factor = area_energy_factors[cfg.EFFECTIVE_AREA_0_30_DEG_INDEX - 1]
    elif dec <= 0 and dec > -5:
        area_energy_factor = area_energy_factors[cfg.EFFECTIVE_AREA_MIN5_0_DEG_INDEX - 1]
    elif dec <= -5 and dec > -30:
        area_energy_factor = area_energy_factors[
            cfg.EFFECTIVE_AREA_MIN30_MIN5_DEG_INDEX - 1
        ]
    elif dec <= -30 and dec >= -90:
        area_energy_factor = area_energy_factors[
            cfg.EFFECTIVE_AREA_MIN90_MIN30_DEG_INDEX - 1
        ]
    constant = (
        (hubble_in_s ** 2) * seconds / (4 * np.pi * (z ** 2) * (cfg.SPEED_OF_LIGHT ** 2))
    )  # m^-2 * s
    expected_nu = constant * nu_flux * (cfg.E0 ** 2) * area_energy_factor
    return expected_nu


def flux_contribute(z, dec, nu_flux=cfg.FLUX_NU):
    """
    Given the redshift and the declination of a source, determines the contribution
    to the test statistic related to the neutrino flux of the source
    """
    mu = expected_nu_from_source(z, dec, nu_flux=nu_flux)
    contribute = np.log(1 - (1 + mu) * np.exp(-mu))
    return contribute


fig, ax = plt.subplots(figsize=(6, 6))

xx = np.logspace(-5, 10, 100)
ress1 = []
ress2 = []
ress3 = []
ress4 = []
for x in xx:
    res1 = flux_contribute(x, cfg.DECS, nu_flux=cfg.FLUX1)
    res2 = flux_contribute(x, cfg.DECS, nu_flux=cfg.FLUX2)
    res3 = flux_contribute(x, cfg.DECS, nu_flux=cfg.FLUX3)
    res4 = flux_contribute(x, cfg.DECS, nu_flux=cfg.FLUX4)
    ress1.append(res1)
    ress2.append(res2)
    ress3.append(res3)
    ress4.append(res4)
plt.plot(
    xx,
    ress4,
    linewidth=cfg.LINEWIDTH_FLUXPLOT,
    label="$\phi_0=10^{%d}$ s$^{-1}$GeV$^{-1}$" % (int(np.log10(cfg.FLUX4))),
)
plt.plot(
    xx,
    ress2,
    linewidth=cfg.LINEWIDTH_FLUXPLOT,
    label="$\phi_0=10^{%d}$ s$^{-1}$GeV$^{-1}$" % (int(np.log10(cfg.FLUX2))),
)
plt.plot(
    xx,
    ress1,
    linewidth=cfg.LINEWIDTH_FLUXPLOT,
    label="$\phi_0=10^{%d}$ s$^{-1}$GeV$^{-1}$" % (int(np.log10(cfg.FLUX1))),
)
plt.xlabel("Redshift", fontsize="x-large")
plt.ylabel("$\log{[1-(1+\mu)e^{-\mu}]}$ [a.u.]", fontsize="x-large")
plt.ylim(-30, 20)

plt.axvline(0.017, color="black", linestyle="--", linewidth=2, label="NGC 7469")
plt.axvline(0.001, color="black", linestyle="dotted", linewidth=2, label="Nearest AGN")
ax.set_xticklabels(ax.get_xticks(), size=14)
ax.set_yticklabels(ax.get_yticks(), size=14)
plt.legend(fontsize="large", loc="upper right")
plt.xscale("log")
plt.title("TS dependency on $\phi_0$", fontsize="xx-large")
plt.xlim(2e-4, 2e2)
plt.savefig(figures_path / "flux_component", bbox_inches="tight")
plt.savefig(figures_path / "flux_component.pdf", bbox_inches="tight")
