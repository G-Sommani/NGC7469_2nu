import numpy as np
import matplotlib.pyplot as plt
from math import factorial
import csv
import config as cfg
from loading_functions import Loader

loader = Loader()
data_path = loader.data_path
figures_path = loader.figures_path

print("Loading NGC7469 SED...")

# Load NGC7469 SED
f_7469_SED = open(data_path / "NGC7469 _SED_with_sources.txt", "r")

energies_7469 = np.array([])
fluxes_7469 = np.array([])
flux_errs_7469 = np.array([])

for x in f_7469_SED:
    if x[0] != "#":
        datas = x.split("   ")
        energy = float(datas[0])
        flux = float(datas[2])
        flux_err = float(datas[3].split(" ")[0])

        energies_7469 = np.append(energies_7469, energy)
        fluxes_7469 = np.append(fluxes_7469, flux)
        flux_errs_7469 = np.append(flux_errs_7469, flux_err)

print("Loading NGC1068 SED...")

# Load NGC1068 SED
f_1068_SED = open(data_path / "NGC1068_SED_with_sources.txt", "r")

energies_1068 = np.array([])
fluxes_1068 = np.array([])
flux_errs_1068 = np.array([])

for x in f_1068_SED:
    if x[0] != "#":
        datas = x.split("   ")
        energy = float(datas[0])
        flux = float(datas[2])
        flux_err = float(datas[3].split(" ")[0])

        energies_1068 = np.append(energies_1068, energy)
        fluxes_1068 = np.append(fluxes_1068, flux)
        flux_errs_1068 = np.append(flux_errs_1068, flux_err)

# Define years of integration
days_alerts_gb = cfg.MJD_04102023 - cfg.MJD_GOLDBRONZE_START
days_alerts_icecat = cfg.MJD_04102023 - cfg.MJD_ICECAT_START
years_alerts_gb = days_alerts_gb / cfg.DAYS_IN_YEAR
years_alerts_icecat = days_alerts_icecat / cfg.DAYS_IN_YEAR

print("Calculating confidence intervals...")


# Confidence intervals calculation
def lr(N, lam, T):
    """
    likelihood ratio for a poissonian likelihood where mu = lam*T
    """
    if N == 0:
        return np.exp(-lam * T)
    y = lam * T / N
    return (y ** N) * np.exp(N - lam * T)


def lr_evaluator(lam, T=years_alerts_gb, until=100):
    """
    Given an event rate 'lam' and a period of time 'T',
    evaluate all the poissonian likelihood ratios for all
    the counts 'N' from 0 until a fixed integer 'until'
    """
    counts_space = np.linspace(1, until, until)
    outcomes = [np.exp(-lam * T)]
    counts_space_final = [int(0)]
    for count in counts_space:
        outcome = lr(count, lam, T)
        outcomes.append(outcome)
        counts_space_final.append(int(count))
    counts_space_final = np.asarray(counts_space_final)
    outcomes = np.asarray(outcomes)
    return counts_space_final, outcomes


def pois(x, lam, T=years_alerts_gb):
    """
    poissonian probability
    """
    x = int(x)
    return (((lam * T) ** x) / factorial(x)) * np.exp(-lam * T)


def find_coverage_interval(lam, CL=0.9, T=years_alerts_gb, until=100):
    """
    Given a confidence level 'CL', it finds, for a specific rate 'lam'
    and time interval 'T', the interval in counts that ensures a
    coverage of at least CL.
    """
    counts_space, outcomes = lr_evaluator(lam, T=T, until=until)
    antithreshold_space = np.linspace(0, 1, 1000)
    probabilities = np.array([])
    inf_count_limit = None
    sup_count_limit = None
    for count in counts_space:
        probability = pois(count, lam, T=T)
        probabilities = np.append(probabilities, probability)
    # This cylce 'for' tries all possible thresholds for the likelihood
    # ratio. To each threshold corresponds only one coverage interval.
    # Once the coverage interval related to the specific threshold
    # is found, it calculates the respective coverage. As soon as
    # the coverage is greater or equal to the confidence level, this
    # is returned.
    for antithreshold in antithreshold_space:
        threshold = 1 - antithreshold
        for index, count in enumerate(counts_space):
            lr_outcome = outcomes[index]
            if lr_outcome >= threshold:
                inf_count_limit = count
                break
        for i in range(len(outcomes)):
            count = counts_space[-1 - i]
            lr_outcome = outcomes[-1 - i]
            if lr_outcome >= threshold:
                sup_count_limit = count
                break
        real_coverage = np.sum(probabilities[outcomes >= threshold])
        if real_coverage >= CL:
            return inf_count_limit, sup_count_limit


# Generate the confidence belt
inf_count_limits_gb = np.array([])
sup_count_limits_gb = np.array([])
inf_count_limits_icecat = np.array([])
sup_count_limits_icecat = np.array([])
rates = np.linspace(0, 200, 201)
rates = rates / 100
for rate in rates:
    inf_count_limit_gb, sup_count_limit_gb = find_coverage_interval(rate, T=years_alerts_gb)
    inf_count_limits_gb = np.append(inf_count_limits_gb, inf_count_limit_gb)
    sup_count_limits_gb = np.append(sup_count_limits_gb, sup_count_limit_gb)
    inf_count_limit_icecat, sup_count_limit_icecat = find_coverage_interval(rate, T=years_alerts_icecat)
    inf_count_limits_icecat = np.append(inf_count_limits_icecat, inf_count_limit_icecat)
    sup_count_limits_icecat = np.append(sup_count_limits_icecat, sup_count_limit_icecat)

print("Plotting confidence belt...")

# Plot confidence belt
plt.fill_between(
    rates, 
    inf_count_limits_gb, 
    sup_count_limits_gb, 
    color="blue", 
    alpha=0.2,
    label="4 yrs",
)
plt.plot(rates, inf_count_limits_gb, color="blue")
plt.plot(rates, sup_count_limits_gb, color="blue")
plt.fill_between(
    rates, 
    inf_count_limits_icecat, 
    sup_count_limits_icecat, 
    color="red", 
    alpha=0.2,
    label = "12 yrs",
)
plt.plot(rates, inf_count_limits_icecat, color="red")
plt.plot(rates, sup_count_limits_icecat, color="red")
plt.xlabel(r"Rate of neutrinos [year$^{-1}$]")
plt.ylabel("Number of detected neutrinos")
plt.title("Confidence belt")
plt.legend()
plt.savefig(figures_path / "confidence_belt")
plt.close()

print(f"Confidence belt plot saved in {figures_path / 'confidence_belt.png'}")

# Given 2 detected neutrinos, the superior limit on the rate
lim_sup_2_nu_gb = max(rates[inf_count_limits_gb == 2.0])
lim_sup_2_nu_icecat = max(rates[inf_count_limits_icecat == 2.0])

# Given 2 detected neutrinos, the inferior limit on the rate
lim_inf_2_nu_gb = min(rates[sup_count_limits_gb == 2.0])
lim_inf_2_nu_icecat = min(rates[sup_count_limits_icecat == 2.0])

print("Loading effective areas...")

# Load effective areas (from IceCat-1: the IceCube Event Catalog of Alert Tracks,
# The IceCube Collaboration, arXiv e-prints, arXiv:2304.01174)
f_effA = open(data_path / "Effa_all_streams_gold_bronze.txt", "r")
# These effective areas are binned in neutrino energy (10 TeV - 500 PeV)
# and in broad declination (stored internally in the files as local
# zenith angle) bins that reflect the major changes in effective area
# due to background rejection in the down-going direction and Earth
# absorption at higher energies for the up-going directions.
# * Declination boundaries (90 deg, 30 deg, 0 deg, -5 deg, -30 deg, -90 deg)
# We will look only in the declination boundary between 30 and 0 degrees
effA2022_found = False
effA2023_found = False
areas = list()
energies = list()
k = 0
for i, x in enumerate(f_effA):
    if i > 0:
        Energy = float(x.split()[0])
        Area = float(x.split()[2])
        energies.append(Energy)
        areas.append(Area)
        k += 1
        # Energy bins that concern the neutrino energies
        if Energy > cfg.E_NU_2022 and not effA2022_found:
            effA2022 = areas[k - 1]
            effA2022_found = True
        if Energy > cfg.E_NU_2023 and not effA2023_found:
            effA2023 = areas[k - 1]
            effA2023_found = True

print("Converting rates to energy fluxes...")

# Convert energy in erg
E22_erg = cfg.E_NU_2022 * cfg.GEV_TO_ERG
E23_erg = cfg.E_NU_2023 * cfg.GEV_TO_ERG

# Convert effective areas in cm^2
A22_cm = effA2022 * 1e4
A23_cm = effA2023 * 1e4

# Average energy over effective area for the two neutrinos
Avg_en_effA = (E22_erg / A22_cm + E23_erg / A23_cm) / cfg.NUMBER_OF_NU

# Convert rate in s^-1
rate_sup_2nu_gb = lim_sup_2_nu_gb / cfg.SECONDS_IN_YEAR
rate_inf_2nu_gb = lim_inf_2_nu_gb / cfg.SECONDS_IN_YEAR
rate_sup_2nu_icecat = lim_sup_2_nu_icecat / cfg.SECONDS_IN_YEAR
rate_inf_2nu_icecat = lim_inf_2_nu_icecat / cfg.SECONDS_IN_YEAR

# Calculate superior limit over energy flux, for the average energy of the two neutrinos (THIS IS NOT PLOTTED)
E_sup_Flux_2_nu_gb = Avg_en_effA * rate_sup_2nu_gb
E_inf_Flux_2_nu_gb = Avg_en_effA * rate_inf_2nu_gb
E_sup_Flux_2_nu_icecat = Avg_en_effA * rate_sup_2nu_icecat
E_inf_Flux_2_nu_icecat = Avg_en_effA * rate_inf_2nu_icecat

# Calculate averaged energy of the two neutrinos
averaged_energy = (cfg.E_NU_2022 * effA2022 + cfg.E_NU_2023 * effA2023) / (
    effA2022 + effA2023
)

# As lower bound in energy we take 100 TeV less than the average energy.
# (Similar to TXS paper)
min_energy = min(cfg.E_NU_2022, cfg.E_NU_2023)
max_energy = max(cfg.E_NU_2022, cfg.E_NU_2023)
lower_energy_bound = min_energy - 1e5

# As superior bound in energy we take 20 PeV (circa the highest energy neutrinos that IceCube can hope to detect).
energy_space_nu_flux_7469_general = np.logspace(
    np.log10(lower_energy_bound), np.log10(cfg.SUPERIOR_ENERGY_BOUND), 500
)
# We also define another energy space just between the two reported neutrino energies.
energy_space_nu_flux_7469_between_nu = np.logspace(
    np.log10(min_energy), np.log10(max_energy), 100
)

# interpolate the effective areas over the various energy spaces and convert them in cm.
effas_interp_general = (
    np.interp(energy_space_nu_flux_7469_general, energies, areas) * cfg.M2_TO_CM2
)
effas_interp_between_nu = (
    np.interp(energy_space_nu_flux_7469_between_nu, energies, areas) * cfg.M2_TO_CM2
)

# estimate superior and lower limits on energies (THESE ARE PLOTTED):
#     1) Over the whole energy uncertainty;
#     2) Between the two reported neutrino energies.
E_sup_Flux_2_nu_gb_general = (
    energy_space_nu_flux_7469_general
    * cfg.GEV_TO_ERG
    * rate_sup_2nu_gb
    / (effas_interp_general)
)
E_inf_Flux_2_nu_gb_general = (
    energy_space_nu_flux_7469_general
    * cfg.GEV_TO_ERG
    * rate_inf_2nu_gb
    / (effas_interp_general)
)
E_sup_Flux_2_nu_gb_between_nu = (
    energy_space_nu_flux_7469_between_nu
    * cfg.GEV_TO_ERG
    * rate_sup_2nu_gb
    / (effas_interp_between_nu)
)
E_inf_Flux_2_nu_gb_between_nu = (
    energy_space_nu_flux_7469_between_nu
    * cfg.GEV_TO_ERG
    * rate_inf_2nu_gb
    / (effas_interp_between_nu)
)

E_sup_Flux_2_nu_icecat_general = (
    energy_space_nu_flux_7469_general
    * cfg.GEV_TO_ERG
    * rate_sup_2nu_icecat
    / (effas_interp_general)
)
E_inf_Flux_2_nu_icecat_general = (
    energy_space_nu_flux_7469_general
    * cfg.GEV_TO_ERG
    * rate_inf_2nu_icecat
    / (effas_interp_general)
)
E_sup_Flux_2_nu_icecat_between_nu = (
    energy_space_nu_flux_7469_between_nu
    * cfg.GEV_TO_ERG
    * rate_sup_2nu_icecat
    / (effas_interp_between_nu)
)
E_inf_Flux_2_nu_icecat_between_nu = (
    energy_space_nu_flux_7469_between_nu
    * cfg.GEV_TO_ERG
    * rate_inf_2nu_icecat
    / (effas_interp_between_nu)
)

print(f"Inferior limit on energy flux 4yrs: {E_inf_Flux_2_nu_gb:.2} erg cm-2 s-1")
print(f"Superior limit on energy flux 4yrs: {E_sup_Flux_2_nu_gb:.2} erg cm-2 s-1")
print(f"\n\n**************************************************\n\n")
print(
    f"90% CL of neutrino flux at {averaged_energy/1e3:.0f} TeV: [{E_inf_Flux_2_nu_gb/(cfg.TEV_TO_ERG*(averaged_energy*1e-3)**2):.3},{E_sup_Flux_2_nu_gb/(cfg.TEV_TO_ERG*(averaged_energy*1e-3)**2):.3}] TeV-1*cm-2*s-1"
)
print(
    f"Lower and superior bounds in energy: {lower_energy_bound/1e3:.0f} TeV, {cfg.SUPERIOR_ENERGY_BOUND/1e6:.0f} PeV"
)
print(f"\n\n**************************************************\n\n")

print(f"Inferior limit on energy flux 12yrs: {E_inf_Flux_2_nu_icecat:.2} erg cm-2 s-1")
print(f"Superior limit on energy flux 12yrs: {E_sup_Flux_2_nu_icecat:.2} erg cm-2 s-1")
print(f"\n\n**************************************************\n\n")
print(
    f"90% CL of neutrino flux at {averaged_energy/1e3:.0f} TeV: [{E_inf_Flux_2_nu_icecat/(cfg.TEV_TO_ERG*(averaged_energy*1e-3)**2):.3},{E_sup_Flux_2_nu_icecat/(cfg.TEV_TO_ERG*(averaged_energy*1e-3)**2):.3}] TeV-1*cm-2*s-1"
)
print(
    f"Lower and superior bounds in energy: {lower_energy_bound/1e3:.0f} TeV, {cfg.SUPERIOR_ENERGY_BOUND/1e6:.0f} PeV"
)
print(f"\n\n**************************************************\n\n")

print("Retrieving neutrino flux intervals for NGC1068...")

# Neutrino flux for NGC1068
energies_1068_nu = np.logspace(
    np.log10(cfg.MIN_E_NU_1068), np.log10(cfg.MAX_E_NU_1068), 100
)


def nu_flux(x, flux, gamma):
    """
    given an energy, get nu flux starting with a flux at 1 TeV
    and a spectral index gamma.
    """
    return flux * (x / 1e3) ** (2 - gamma)


nu_fluxes_min_1068 = nu_flux(energies_1068_nu, cfg.FLUX_NU_MIN_1068, cfg.GAMMA_MAX_1068)
nu_fluxes_max_1068 = nu_flux(energies_1068_nu, cfg.FLUX_NU_MAX_1068, cfg.GAMMA_MIN_1068)

print("Retrievieng IceCube 7 yrs sensitivities...")

# Load IceCube 7 yrs sensitivities. Data from
# "All-sky Search for Time-integrated Neutrino
# Emission from Astrophysical Sources with 7 yr
# of IceCube Data", The IceCube Collaboration,
# 2017, AJ, 835, 151. Precisely, from Figure 5.
with open(data_path / "pstracks_7yrs_diffsens.csv", "r") as file:
    csvreader = csv.reader(file)
    fluxes_sens = np.array([])
    for row in csvreader:
        flux_TeV = float((row[0] + "." + row[1]).split()[1])  # Unit is TeV*cm-2*s-1
        flux = flux_TeV * 1.60218  # Unit: erg*cm-2*s-1
        fluxes_sens = np.append(fluxes_sens, flux)
max_E_7yrs = cfg.MAX_E_7YRS * 1e3  # GeV
min_E_7yrs = cfg.MIN_E_7YRS * 1e3  # GeV
first_bin = min_E_7yrs * 10 ** (cfg.LOG_WIDTH_7YRS_BINS / 2)
fluxes_7yrs = [fluxes_sens[0]]
bins_7yrs = [min_E_7yrs]
for i in range(14):
    new_bin = bins_7yrs[-1] * 10 ** (cfg.LOG_WIDTH_7YRS_BINS)
    if i != 13:
        bins_7yrs.append(new_bin)
        bins_7yrs.append(new_bin)
        fluxes_7yrs.append(fluxes_sens[i])
        fluxes_7yrs.append(fluxes_sens[i + 1])
    else:
        bins_7yrs.append(new_bin)
        fluxes_7yrs.append(fluxes_sens[i])

print("Retrieving intrinsic flux data...")


# Retrieve intrinsic fluxes. Data from "BAT AGN Spectroscopic Survey.
# V. X-Ray Properties of the Swift/BAT 70-month AGN Catalog", C. Ricci
# et al., 2017, ApJS 233 17.
def find_fluxes(
    gamma,  # Spectral index
    E2phi_tot,  # Flux
    gamma_err,  # Spectral index err
    E2phi_tot_err,  # Flux err
    E_1=30,  # Energy first point in SED
    E_2=70,  # Energy second point in SED
    E_min=14,  # Minimum energy represented by flux+spectral index
    E_MAX=195,  # Maximum energy represented by flux+spectral index
):
    """
    Function to translate a flux + spectral index into two fluxes to insert
    in an SED at the energies E_1 and E_2. Here specifically used to retrieve
    hard x-ray data.
    """
    phi_1 = ((2 - gamma) * (E_1 ** (-gamma)) * E2phi_tot) / (
        E_MAX ** (2 - gamma) - E_min ** (2 - gamma)
    )
    phi_2 = phi_1 * ((E_1 / E_2) ** gamma)
    dphi_1dgamma = -phi_1 * (
        1 / (2 - gamma)
        + np.log(E_1)
        - (
            np.log(E_MAX) * (E_MAX ** (2 - gamma))
            - np.log(E_min) * (E_min ** (2 - gamma))
        )
        / (E_MAX ** (2 - gamma) - E_min ** (2 - gamma))
    )
    dphi_1dphi_tot = phi_1 / E2phi_tot
    dphi_2dgamma = phi_2 * np.log(E_1 / E_2)
    dphi_2dphi_1 = phi_2 / phi_1
    phi_1_err = np.sqrt(
        (gamma_err * dphi_1dgamma) ** 2 + (E2phi_tot_err * dphi_1dphi_tot) ** 2
    )
    phi_2_err = np.sqrt(
        (gamma_err * dphi_2dgamma) ** 2 + (phi_1_err * dphi_2dphi_1) ** 2
    )
    E2phi_1 = phi_1 * (E_1 ** 2)
    E2phi_2 = phi_2 * (E_2 ** 2)
    E2phi_1_err = phi_1_err * (E_1 ** 2)
    E2phi_2_err = phi_2_err * (E_2 ** 2)
    return E2phi_1, E2phi_2, E2phi_1_err, E2phi_2_err


intrinsic_data_1068 = find_fluxes(
    cfg.INTRINSIC_GAMMA_1068,
    cfg.INTRINSIC_HARD_XRAY_FLUX_1068,
    cfg.INTRINSIC_GAMMA_ERR_1068,
    cfg.INTRINSIC_HARD_XRAY_FLUX_ERR_1068,
)
intrinsic_data_7469 = find_fluxes(
    cfg.INTRINSIC_GAMMA_7469,
    cfg.INTRINSIC_HARD_XRAY_FLUX_7469,
    cfg.INTRINSIC_GAMMA_ERR_7469,
    cfg.INTRINSIC_HARD_XRAY_FLUX_ERR_7469,
)
intr_fluxes_sed_1068 = intrinsic_data_1068[:2]
intr_fluxes_sed_errs_1068 = intrinsic_data_1068[2:]
intr_fluxes_sed_7469 = intrinsic_data_7469[:2]
intr_fluxes_sed_errs_7469 = intrinsic_data_7469[2:]

print("Start plotting SED...")

# SED plot
fig, ax = plt.subplots(figsize=(9, 4.5))
plt.errorbar(
    energies_7469,
    fluxes_7469,
    flux_errs_7469,
    linestyle="",
    color="black",
    alpha=0.8,
    marker=".",
    markersize=3,
    label="NGC 7469, electromagnetic observations",
)
plt.errorbar(
    energies_1068,
    fluxes_1068,
    flux_errs_1068,
    linestyle="",
    alpha=0.8,
    marker=".",
    markersize=3,
    label="NGC 1068, electromagnetic observations",
    color="darkgrey",
)
plt.fill_between(
    energy_space_nu_flux_7469_general,
    E_sup_Flux_2_nu_gb_general,
    E_inf_Flux_2_nu_gb_general,
    alpha=0.3,
    color="white",
    hatch="////",
    edgecolor="tab:blue",
    label=r"NGC 7469 $\nu_\mu+\bar{\nu}_\mu$, 90% CL, this work, energy uncertainty",
)
plt.fill_between(
    energy_space_nu_flux_7469_between_nu,
    E_sup_Flux_2_nu_gb_between_nu,
    E_inf_Flux_2_nu_gb_between_nu,
    alpha=1,
    color="darkblue",
    label=r"NGC 7469 $\nu_\mu+\bar{\nu}_\mu$, 90% CL, this work, between neutrino energies",
)
plt.fill_between(
    energies_1068_nu,
    nu_fluxes_min_1068,
    nu_fluxes_max_1068,
    alpha=1,
    label=r"NGC 1068 $\nu_\mu+\bar{\nu}_\mu$, IceCube 2022",
    color="lightgrey",
)
plt.plot(
    bins_7yrs,
    fluxes_7yrs,
    label="Differential sensitivities, IceCube 7 years",
    color="orange",
    linestyle="--",
)
plt.errorbar(
    [cfg.ENERGY1_INTRINSIC_FLUX * 1e-6, cfg.ENERGY2_INTRINSIC_FLUX * 1e-6],
    intr_fluxes_sed_7469,
    intr_fluxes_sed_errs_7469,
    markersize=1,
    linewidth=1,
    label="NGC 7469, hard x-ray (14 - 195 keV) intrinsic flux",
)
plt.errorbar(
    [cfg.ENERGY1_INTRINSIC_FLUX * 1e-6, cfg.ENERGY2_INTRINSIC_FLUX * 1e-6],
    intr_fluxes_sed_1068,
    intr_fluxes_sed_errs_1068,
    markersize=1,
    linewidth=1,
    label="NGC 1068, hard x-ray (14 - 195 keV) intrinsic flux",
)
plt.yscale("log")
plt.xscale("log")
plt.ylim(1e-15, 3e-9)
plt.legend(fontsize="small")
plt.xlabel("Energy [GeV]")
plt.ylabel(r"$E^2\Phi$ [erg cm$^{-2}$ s$^{-1}$]")
plt.savefig(figures_path / "SED_nu_flux", bbox_inches="tight", dpi=200)
plt.savefig(figures_path / "SED_nu_flux.pdf", bbox_inches="tight", dpi=200)
plt.close()

print(f"SED saved in {figures_path / 'SED_nu_flux.png'}.")

print(f"Plotting the estimated flux for 12 yrs...")

plt.subplots(figsize=cfg.FIGSIZE_TS)
plt.fill_between(
    energy_space_nu_flux_7469_general,
    E_sup_Flux_2_nu_gb_general,
    E_inf_Flux_2_nu_gb_general,
    alpha=0.3,
    color="white",
    hatch="////",
    edgecolor="tab:blue",
    label=r"4 yrs, energy uncertainty",
)
plt.fill_between(
    energy_space_nu_flux_7469_between_nu,
    E_sup_Flux_2_nu_gb_between_nu,
    E_inf_Flux_2_nu_gb_between_nu,
    alpha=0.9,
    color="darkblue",
    label=r"4 yrs, between neutrino energies",
)
plt.fill_between(
    energy_space_nu_flux_7469_general,
    E_sup_Flux_2_nu_icecat_general,
    E_inf_Flux_2_nu_icecat_general,
    alpha=0.3,
    color="white",
    hatch="\\\\\\\\",
    edgecolor="tab:red",
    label=r"12 yrs, energy uncertainty",
)
plt.fill_between(
    energy_space_nu_flux_7469_between_nu,
    E_sup_Flux_2_nu_icecat_between_nu,
    E_inf_Flux_2_nu_icecat_between_nu,
    alpha=0.9,
    color="darkred",
    label=r"12 yrs, between neutrino energies",
)

plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.xlabel("Energy [GeV]")
plt.ylabel(r"$E^2\Phi$ [erg cm$^{-2}$ s$^{-1}$]")
plt.title(r"NGC 7469 $\nu_\mu+\bar{\nu}_\mu$, 90% CL")
plt.savefig(figures_path / "nu_flux_12yrs", bbox_inches="tight", dpi=200)
plt.savefig(figures_path / "nu_flux_12yrs.pdf", bbox_inches="tight", dpi=200)
plt.close()