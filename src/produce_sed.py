import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os
from math import factorial
import csv

# Definition of constants
MJD_04102023 = 60132
MJD_GOLDBRONZE_START = 58635
DAYS_IN_YEAR = 365
E_NU_2022 = 1.8399e05  # GeV
E_NU_2023 = 1.2729e05  # GeV
SUPERIOR_ENERGY_BOUND = 2e7 # GeV
GEV_TO_ERG = 1.60218e-3
TEV_TO_ERG = 1e3 * GEV_TO_ERG
SECONDS_IN_YEAR = DAYS_IN_YEAR * 24 * 60 * 60
NUMBER_OF_NU = 2
FLUX_NU_MIN_1068 = 2.9e-11 * 1.6  # erg*cm-2*s-1, at 1 TeV
FLUX_NU_MAX_1068 = 7.1e-11 * 1.6  # erg*cm-2*s-1, at 1 TeV
GAMMA_MIN_1068 = 3.0  # spectral index
GAMMA_MAX_1068 = 3.4  # spectral index
MIN_E_NU_1068 = 1.5e3  # GeV
MAX_E_NU_1068 = 1.5e4  # GeV
MAX_E_7YRS = 1e6  # TeV
MIN_E_7YRS = 1e-1  # TeV
LOG_WIDTH_7YRS_BINS = 0.5
INTRINSIC_HARD_XRAY_FLUX_1068 = 206.00e-12  # erg*cm-2*s-1
INTRINSIC_HARD_XRAY_FLUX_ERR_1068 = 22e-12
INTRINSIC_GAMMA_1068 = 2.37
INTRINSIC_GAMMA_ERR_1068 = 0.1
INTRINSIC_HARD_XRAY_FLUX_7469 = 64.30e-12  # erg*cm-2*s-1
INTRINSIC_HARD_XRAY_FLUX_ERR_7469 = 8e-12
INTRINSIC_GAMMA_7469 = 2.09
INTRINSIC_GAMMA_ERR_7469 = 0.12
ENERGY1_INTRINSIC_FLUX = 30  # keV
ENERGY2_INTRINSIC_FLUX = 70  # keV

print("Definition of paths...")

# Definition of paths
cwd = Path(os.getcwd())
data_path = cwd / "../data"
figures_path = cwd / "../figures"

print(f"Current dir: {cwd}")
print(f"Data dir: {data_path}")
print(f"Figures dir: {figures_path}")

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
days_alerts = MJD_04102023 - MJD_GOLDBRONZE_START
years_alerts = days_alerts / DAYS_IN_YEAR

print("Calculating confidence intervals...")


# Confidence intervals calculation
def lr(N, lam, T):
    """
    likelihood ratio for a poissonian likelihood where mu = lam*T
    """
    if N == 0:
        return np.exp(-lam * T)
    y = lam * T / N
    return (y**N) * np.exp(N - lam * T)


def lr_evaluator(lam, T=years_alerts, until=100):
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


def pois(x, lam, T=years_alerts):
    """
    poissonian probability
    """
    x = int(x)
    return (((lam * T) ** x) / factorial(x)) * np.exp(-lam * T)


def find_coverage_interval(lam, CL=0.9, T=years_alerts, until=100):
    """
    Given a confidence level 'CL', it finds, for a specific rate 'lam'
    and time interval 'T', the interval in counts that ensures a
    coverage of at least CL.
    """
    counts_space, outcomes = lr_evaluator(lam, T=T, until=until)
    antithreshold_space = np.linspace(0, 1, 1000)
    probabilities = np.array([])
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
inf_count_limits = np.array([])
sup_count_limits = np.array([])
rates = np.linspace(0, 200, 201)
rates = rates / 100
for rate in rates:
    inf_count_limit, sup_count_limit = find_coverage_interval(rate, T=years_alerts)
    inf_count_limits = np.append(inf_count_limits, inf_count_limit)
    sup_count_limits = np.append(sup_count_limits, sup_count_limit)

print("Plotting confidence belt...")

# Plot confidence belt
plt.fill_between(rates, inf_count_limits, sup_count_limits, color="blue", alpha=0.2)
plt.plot(rates, inf_count_limits, color="blue")
plt.plot(rates, sup_count_limits, color="blue")
plt.xlabel(r"Rate of neutrinos [year$^{-1}$]")
plt.ylabel("Number of detected neutrinos")
plt.title("Confidence belt")
plt.savefig(figures_path / "confidence_belt")
plt.close()

print(f"Confidence belt plot saved in {figures_path / 'confidence_belt.png'}")

# Given 2 detected neutrinos, the superior limit on the rate
lim_sup_2_nu = max(rates[inf_count_limits == 2.0])

# Given 2 detected neutrinos, the inferior limit on the rate
lim_inf_2_nu = min(rates[sup_count_limits == 2.0])

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
        if Energy > E_NU_2022 and not effA2022_found:
            effA2022 = areas[k - 2]
            effA2022_found = True
        if Energy > E_NU_2023 and not effA2023_found:
            effA2023 = areas[k - 2]
            effA2023_found = True

print("Converting rates to energy fluxes...")

# Convert energy in erg
E22_erg = E_NU_2022 * GEV_TO_ERG
E23_erg = E_NU_2023 * GEV_TO_ERG

# Convert effective areas in cm^2
A22_cm = effA2022 * 1e4
A23_cm = effA2023 * 1e4

# Average energy over effective area for the two neutrinos
Avg_en_effA = (E22_erg / A22_cm + E23_erg / A23_cm) / NUMBER_OF_NU

# Convert rate in s^-1
rate_sup_2nu = lim_sup_2_nu / SECONDS_IN_YEAR
rate_inf_2nu = lim_inf_2_nu / SECONDS_IN_YEAR

# Calculate superior limit over energy flux
E_sup_Flux_2_nu = Avg_en_effA * rate_sup_2nu
E_inf_Flux_2_nu = Avg_en_effA * rate_inf_2nu

# Translate limits in something easy to plot
half_log_space_2_nu = 10 ** (np.log10(E_inf_Flux_2_nu * E_sup_Flux_2_nu) / 2)
yerr_2_nu = [
    [half_log_space_2_nu - E_inf_Flux_2_nu],
    [E_sup_Flux_2_nu - half_log_space_2_nu],
]

# Calculate averaged energy of the two neutrinos
averaged_energy = (E_NU_2022 * effA2022 + E_NU_2023 * effA2023) / (effA2022 + effA2023)

# As lower bound in energy we take 100 TeV less than the average energy.
# (Similar to TXS paper)
min_energy = min(E_NU_2022, E_NU_2023)
lower_energy_bound = min_energy - 1e5

# As superior bound in energy we take 20 PeV (circa the highest energy neutrinos that IceCube can hope to detect)
energy_space_nu_flux_7469 = np.logspace(
    1, 12, 100
)
#energy_space_nu_flux_7469 = np.logspace(
#    np.log10(lower_energy_bound), np.log10(SUPERIOR_ENERGY_BOUND), 100
#)

effas_interp = np.interp(energy_space_nu_flux_7469, energies, areas)*1e4

#effas = np.logspace(np.log10(min(A22_cm, A23_cm)), np.log10(areas[-4]*1e4), 100)
E_sup_Flux_2_nu_alt = energy_space_nu_flux_7469*GEV_TO_ERG*rate_sup_2nu/(effas_interp)
E_inf_Flux_2_nu_alt = energy_space_nu_flux_7469*GEV_TO_ERG*rate_inf_2nu/(effas_interp)

print(
    f"Inferior limit on energy flux: {E_inf_Flux_2_nu:.2} erg cm-2 s-1, {E_inf_Flux_2_nu:.2}"
)
print(f"Superior limit on energy flux: {E_sup_Flux_2_nu:.2} erg cm-2 s-1")
print(f"\n\n**************************************************\n\n")
print(
    f"90% CL of neutrino flux at {averaged_energy/1e3:.0f} TeV: [{E_inf_Flux_2_nu/(TEV_TO_ERG*(averaged_energy*1e-3)**2):.3},{E_sup_Flux_2_nu/(TEV_TO_ERG*(averaged_energy*1e-3)**2):.3}] TeV-1*cm-2*s-1"
)
print(
    f"Lower and superior bounds in energy: {lower_energy_bound/1e3:.0f} TeV, {SUPERIOR_ENERGY_BOUND/1e6:.0f} PeV"
)
print(f"\n\n**************************************************\n\n")

print("Retrieving neutrino flux intervals for NGC1068...")

# Neutrino flux for NGC1068
energies_1068_nu = np.logspace(np.log10(MIN_E_NU_1068), np.log10(MAX_E_NU_1068), 100)


def nu_flux(x, flux, gamma):
    """
    given an energy, get nu flux starting with a flux at 1 TeV
    and a spectral index gamma.
    """
    return flux * (x / 1e3) ** (2 - gamma)


nu_fluxes_min_1068 = nu_flux(energies_1068_nu, FLUX_NU_MIN_1068, GAMMA_MAX_1068)
nu_fluxes_max_1068 = nu_flux(energies_1068_nu, FLUX_NU_MAX_1068, GAMMA_MIN_1068)

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
max_E_7yrs = MAX_E_7YRS * 1e3  # GeV
min_E_7yrs = MIN_E_7YRS * 1e3  # GeV
first_bin = min_E_7yrs * 10 ** (LOG_WIDTH_7YRS_BINS / 2)
fluxes_7yrs = [fluxes_sens[0]]
bins_7yrs = [min_E_7yrs]
for i in range(14):
    new_bin = bins_7yrs[-1] * 10 ** (LOG_WIDTH_7YRS_BINS)
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
    E2phi_1 = phi_1 * (E_1**2)
    E2phi_2 = phi_2 * (E_2**2)
    E2phi_1_err = phi_1_err * (E_1**2)
    E2phi_2_err = phi_2_err * (E_2**2)
    return E2phi_1, E2phi_2, E2phi_1_err, E2phi_2_err


intrinsic_data_1068 = find_fluxes(
    INTRINSIC_GAMMA_1068,
    INTRINSIC_HARD_XRAY_FLUX_1068,
    INTRINSIC_GAMMA_ERR_1068,
    INTRINSIC_HARD_XRAY_FLUX_ERR_1068,
)
intrinsic_data_7469 = find_fluxes(
    INTRINSIC_GAMMA_7469,
    INTRINSIC_HARD_XRAY_FLUX_7469,
    INTRINSIC_GAMMA_ERR_7469,
    INTRINSIC_HARD_XRAY_FLUX_ERR_7469,
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
"""
plt.fill_between(
    energy_space_nu_flux_7469,
    half_log_space_2_nu - yerr_2_nu[0][0],
    half_log_space_2_nu + yerr_2_nu[1][0],
    alpha=0.3,
    label=r"NGC 7469 $\nu_\mu+\bar{\nu}_\mu$, 90% CL, this work",
)
"""
plt.fill_between(
    energy_space_nu_flux_7469,
    E_sup_Flux_2_nu_alt,
    E_inf_Flux_2_nu_alt,
    alpha=0.3,
    label=r"NGC 7469 $\nu_\mu+\bar{\nu}_\mu$, 90% CL, this work",
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
    [ENERGY1_INTRINSIC_FLUX * 1e-6, ENERGY2_INTRINSIC_FLUX * 1e-6],
    intr_fluxes_sed_7469,
    intr_fluxes_sed_errs_7469,
    markersize=1,
    linewidth=1,
    label="NGC 7469, hard x-ray (14 - 195 keV) intrinsic flux",
)
plt.errorbar(
    [ENERGY1_INTRINSIC_FLUX * 1e-6, ENERGY2_INTRINSIC_FLUX * 1e-6],
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
plt.close()

print(f"SED saved in {figures_path / 'SED_nu_flux.png'}.")
