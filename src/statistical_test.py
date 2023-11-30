import sys
from pathlib import Path
import os
import requests
import zipfile
import pandas as pd
import numpy as np

ALLOWED_CATALOGS = ["turin", "milliquas"]
ALLOWED_RECONSTRUCTIONS = ["splinempe", "millipede"]
TURIN_INDEX = 0
MILLIQUAS_INDEX = 1
SPLINEMPE_INDEX = 0
MILLIPEDE_INDEX = 1
DEFAULT_CATALOG = ALLOWED_CATALOGS[TURIN_INDEX]
DEFAULT_RECO = ALLOWED_RECONSTRUCTIONS[SPLINEMPE_INDEX]
MILLIQUAS_FILENAME = "milliquas.txt"
MILLIQUAS_URL = "https://quasars.org/milliquas.zip"
MILLIQUAS_ZIP = "milliquas.zip"
MILLIQUAS_RA = "RA"
MILLIQUAS_DEC = "DEC"
MILLIQUAS_Z = "Z"
MILLIQUAS_TYPE = "Type"
MILLIQUAS_COLSPECS = [
    (0, 11),
    (12, 23),
    (25, 50),
    (51, 55),
    (56, 61),
    (62, 67),
    (68, 70),
    (71, 72),
    (73, 74),
    (76, 82),
    (83, 88),
    (90, 95),
    (97, 99),
    (101, 103),
    (105, 126),
    (128, 149),
    (151, 172),
    (174, 195),
]
MILLIQUAS_HEADER = [
    MILLIQUAS_RA,
    MILLIQUAS_DEC,
    "Name",
    MILLIQUAS_TYPE,
    "Rmag",
    "Bmag",
    "Comment",
    "R",
    "B",
    MILLIQUAS_Z,
    "Cite",
    "Zcite",
    "RXpct",
    "Qpct",
    "Xname",
    "Rname",
    "Lobe1",
    "Lobe2",
]
MILLIQUAS_Z_CUT = 2
MILLIQUAS_AGN_CATEGORIES = [
    "ARX",
    "AX",
    "QRX",
    "QX",
    "BRX",
    "BX",
    "KRX",
    "KX",
    "NRX",
    "NX",
]
TURIN_FILENAME = "Turin_Catalogue_Table2_SourceProperties.txt"
TURIN_URL = "https://cdsarc.cds.unistra.fr/ftp/J/A+A/659/A32/tablea2.dat"
TURIN_RAh = "RAh"
TURIN_RAm = "RAm"
TURIN_RAs = "RAs"
TURIN_DEd = "DEd"
TURIN_DEm = "DEm"
TURIN_DEs = "DEs"
TURIN_z = "z"
TURIN_HEADER = [
    "ID",
    "SYCAT",
    "WISE",
    TURIN_RAh,
    TURIN_RAm,
    TURIN_RAs,
    TURIN_DEd,
    TURIN_DEm,
    TURIN_DEs,
    TURIN_z,
    "LD",
    "Class",
    "r_z",
    "ClassL",
    "r_ClassL1",
    "r_ClassL2",
]
TURIN_HEADER_PRESENT = None
EFFECTIVE_AREA_FILENAME = "Effa_all_streams_gold_bronze.txt"
EFFECTIVE_AREA_ENERGY_BINS_INDEX = 0
EFFECTIVE_AREA_30_90_DEG_INDEX = 1
EFFECTIVE_AREA_0_30_DEG_INDEX = 2
EFFECTIVE_AREA_MIN5_0_DEG_INDEX = 3
EFFECTIVE_AREA_MIN30_MIN5_DEG_INDEX = 4
EFFECTIVE_AREA_MIN90_MIN30_DEG_INDEX = 5
HUBBLE_CONSTANT = 73  # km s^-1 Mpc^-1
MPC_TO_KM = 3.086e19
SPEED_OF_LIGHT = 299792458  # m s^-1
MJD_07102023 = 60135
MJD_GOLDBRONZE_START = 58635
DAYS_IN_YEAR = 365
SECONDS_IN_YEAR = DAYS_IN_YEAR * 24 * 60 * 60
FLUX_NU = 1e36  # s-1 GeV-1 at 100 GeV
E0 = 100  # GeV

input_list = sys.argv[1:]


def check_input_name(string: str, catalog: str, reco: str):
    if string in ALLOWED_CATALOGS:
        catalog = string
    elif string in ALLOWED_RECONSTRUCTIONS:
        reco = string
    else:
        print(f"\n******************\n\nAllowed catalogs: {ALLOWED_CATALOGS}")
        print(
            f"Allowed reconstructions: {ALLOWED_RECONSTRUCTIONS}\n\n******************\n"
        )
        raise NameError(f"'{string}' is an incorrect input")
    return catalog, reco


catalog = DEFAULT_CATALOG
reco = DEFAULT_RECO
if len(input_list) == 1:
    input_string = input_list[0]
    catalog, reco = check_input_name(input_string, catalog, reco)
elif len(input_list) == 2:
    input_1 = input_list[0]
    input_2 = input_list[1]
    if input_1 in ALLOWED_CATALOGS and input_2 in ALLOWED_CATALOGS:
        raise Exception(f"Possible to specify the catalog only once")
    if input_1 in ALLOWED_RECONSTRUCTIONS and input_2 in ALLOWED_RECONSTRUCTIONS:
        raise Exception(f"Possible to specify the reconstruction only once")
    catalog, reco = check_input_name(input_1, catalog, reco)
    catalog, reco = check_input_name(input_2, catalog, reco)
elif len(input_list) > 2:
    raise Exception(f"Too many inputs")

print("Definition of paths...")

# Definition of paths
cwd = Path(os.getcwd())
data_path = cwd / "../data"
figures_path = cwd / "../figures"

if catalog == ALLOWED_CATALOGS[TURIN_INDEX]:
    filename = TURIN_FILENAME
    url = TURIN_URL
elif catalog == ALLOWED_CATALOGS[MILLIQUAS_INDEX]:
    filename = MILLIQUAS_FILENAME
    url = MILLIQUAS_URL

print(f"Checking if '{filename}' is in '{data_path}'...")

if os.path.isfile(data_path / filename):
    print(f"'{filename}' in '{data_path}', no need to download")

else:
    print(f"{filename} not found, download {catalog} catalog...")

    r = requests.get(url, allow_redirects=True)
    if catalog == ALLOWED_CATALOGS[TURIN_INDEX]:
        catalog_path = data_path / TURIN_FILENAME
        open(catalog_path, "wb").write(r.content)
    elif catalog == ALLOWED_CATALOGS[MILLIQUAS_INDEX]:
        zip_path = data_path / MILLIQUAS_ZIP
        open(zip_path, "wb").write(r.content)

        print(f"Unzipping {catalog} catalog...")

        with zipfile.ZipFile(zip_path, "r") as zip_ref:
            zip_ref.extractall(data_path)

        print(f"Removing {MILLIQUAS_ZIP}...")

        os.remove(zip_path)

print(f"Loading the {catalog} catalog...")

if catalog == ALLOWED_CATALOGS[TURIN_INDEX]:
    dataframe = pd.read_fwf(
        data_path / TURIN_FILENAME, header=TURIN_HEADER_PRESENT, names=TURIN_HEADER
    )

    def hms_to_deg(h, m, s):
        return 15.0 * (h + (m + s / 60.0) / 60.0)

    def dms_to_deg(d, m, s):
        return d + (m + s / 60.0) / 60.0

    RAs_catalog = hms_to_deg(
        dataframe[TURIN_RAh], dataframe[TURIN_RAm], dataframe[TURIN_RAs]
    ).to_numpy()
    DECs_catalog = dms_to_deg(
        dataframe[TURIN_DEd], dataframe[TURIN_DEm], dataframe[TURIN_DEs]
    ).to_numpy()
    redshifts_catalog = dataframe[TURIN_z].to_numpy()
elif catalog == ALLOWED_CATALOGS[MILLIQUAS_INDEX]:
    dataframe = pd.read_fwf(
        data_path / MILLIQUAS_FILENAME,
        colspecs=MILLIQUAS_COLSPECS,
        names=MILLIQUAS_HEADER,
    )

    print(f"Selecting only AGN within a redshift of {MILLIQUAS_Z_CUT}..")

    # Consider only the agn in the milliquas catalog
    mask01 = np.ma.mask_or(
        dataframe[MILLIQUAS_TYPE] == MILLIQUAS_AGN_CATEGORIES[0],
        dataframe[MILLIQUAS_TYPE] == MILLIQUAS_AGN_CATEGORIES[1],
    )
    mask23 = np.ma.mask_or(
        dataframe[MILLIQUAS_TYPE] == MILLIQUAS_AGN_CATEGORIES[2],
        dataframe[MILLIQUAS_TYPE] == MILLIQUAS_AGN_CATEGORIES[3],
    )
    mask45 = np.ma.mask_or(
        dataframe[MILLIQUAS_TYPE] == MILLIQUAS_AGN_CATEGORIES[4],
        dataframe[MILLIQUAS_TYPE] == MILLIQUAS_AGN_CATEGORIES[5],
    )
    mask67 = np.ma.mask_or(
        dataframe[MILLIQUAS_TYPE] == MILLIQUAS_AGN_CATEGORIES[6],
        dataframe[MILLIQUAS_TYPE] == MILLIQUAS_AGN_CATEGORIES[7],
    )
    mask89 = np.ma.mask_or(
        dataframe[MILLIQUAS_TYPE] == MILLIQUAS_AGN_CATEGORIES[8],
        dataframe[MILLIQUAS_TYPE] == MILLIQUAS_AGN_CATEGORIES[9],
    )
    mask0123 = np.ma.mask_or(mask01, mask23)
    mask4567 = np.ma.mask_or(mask45, mask67)
    mask01234567 = np.ma.mask_or(mask0123, mask4567)
    mask_agn = np.ma.mask_or(mask01234567, mask89)
    dataframe_agn = dataframe[mask_agn]
    # Sources that are too far away are not significant for the test and slow down the program
    dataframe_nearby = dataframe_agn[dataframe_agn[MILLIQUAS_Z] < MILLIQUAS_Z_CUT]
    # There are two blazars which are reported with a redshift of zero. This is unphyisical and
    # these two sources are therefore removed.
    dataframe_nearby_no0 = dataframe_nearby[dataframe_nearby[MILLIQUAS_Z] != 0]
    RAs_catalog = dataframe_nearby_no0[MILLIQUAS_RA].to_numpy()
    DECs_catalog = dataframe_nearby_no0[MILLIQUAS_DEC].to_numpy()
    redshifts_catalog = dataframe_nearby_no0[MILLIQUAS_Z].to_numpy()

print("Defining the test statistic...")

effective_area = np.genfromtxt(data_path / EFFECTIVE_AREA_FILENAME)
energy_bins = effective_area[:, EFFECTIVE_AREA_ENERGY_BINS_INDEX]
effective_area_array = np.array(
    [
        effective_area[:, EFFECTIVE_AREA_30_90_DEG_INDEX],
        effective_area[:, EFFECTIVE_AREA_0_30_DEG_INDEX],
        effective_area[:, EFFECTIVE_AREA_MIN5_0_DEG_INDEX],
        effective_area[:, EFFECTIVE_AREA_MIN30_MIN5_DEG_INDEX],
        effective_area[:, EFFECTIVE_AREA_MIN90_MIN30_DEG_INDEX],
    ]
)


def energy_factor(bin_index):
    """
    Estimate energy factor for expected number of detected neutrinos
    """
    E_max = energy_bins[bin_index + 1]
    E_min = energy_bins[bin_index]
    factor = (E_max - E_min) / (E_max * E_min)
    return factor


def area_energy_factor(A_index):
    """
    Estimate area and energy factor for expected number of detected neutrinos
    """
    effective_area = effective_area_array[A_index]
    factor = 0
    for i in range(len(energy_bins) - 1):
        element = effective_area[i] * energy_factor(i)
        factor += element
    return factor  # units: m^2 GeV^-1


area_energy_factors = np.array([])
for i in range(len(effective_area_array)):
    area_energy_factors = np.append(area_energy_factors, area_energy_factor(i))
hubble_in_s = HUBBLE_CONSTANT / MPC_TO_KM
days = MJD_07102023 - MJD_GOLDBRONZE_START
seconds = (days / DAYS_IN_YEAR) * SECONDS_IN_YEAR


def expected_nu_from_source(z, dec):
    """
    Given the redshift and the declination of a source, determines the total
    number of expected neutrinos from the source
    """
    if dec <= 90 and dec > 30:
        area_energy_factor = area_energy_factors[EFFECTIVE_AREA_30_90_DEG_INDEX - 1]
    elif dec <= 30 and dec > 0:
        area_energy_factor = area_energy_factors[EFFECTIVE_AREA_0_30_DEG_INDEX - 1]
    elif dec <= 0 and dec > -5:
        area_energy_factor = area_energy_factors[EFFECTIVE_AREA_MIN5_0_DEG_INDEX - 1]
    elif dec <= -5 and dec > -30:
        area_energy_factor = area_energy_factors[
            EFFECTIVE_AREA_MIN30_MIN5_DEG_INDEX - 1
        ]
    elif dec <= -30 and dec >= -90:
        area_energy_factor = area_energy_factors[
            EFFECTIVE_AREA_MIN90_MIN30_DEG_INDEX - 1
        ]
    constant = (
        (hubble_in_s**2) * seconds / (4 * np.pi * (z**2) * (SPEED_OF_LIGHT**2))
    )  # m^-2 * s
    expected_nu = constant * FLUX_NU * (E0**2) * area_energy_factor
    return expected_nu


def flux_contribute(z, dec):
    """
    Given the redshift and the declination of a source, determines the contribution
    to the test statistic related to the neutrino flux of the source
    """
    mu = expected_nu_from_source(z, dec)
    contribute = (
        np.log(0.5) + 2 * np.log(mu) - mu
    )  # Here we assume the limit of low fluxes as valid
    return contribute


def unc_contribute(sigma1, sigma2):
    """
    Contribute to the test statistic related only to the uncertainties of the alerts
    """
    return -2 * np.log(sigma1 * sigma2)


print(flux_contribute(0.001, -2))
