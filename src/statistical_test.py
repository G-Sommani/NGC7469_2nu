import sys
from pathlib import Path
import os
import requests
import zipfile
import pandas as pd
import numpy as np
import time

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
MILLIQUAS_NAME = "Name"
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
    MILLIQUAS_NAME,
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
MILLIQUAS_Z_CUT = 0.5
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
TURIN_NAMES_FILENAME = "Turin_Catalogue_Table1_SourceNames.txt"
TURIN_URL = "https://cdsarc.cds.unistra.fr/ftp/J/A+A/659/A32/tablea2.dat"
TURIN_NAMES_URL = "https://cdsarc.cds.unistra.fr/ftp/J/A+A/659/A32/tablea1.dat"
TURIN_ID = "ID"
TURIN_SWIFT = "BAT105-Swift"
TURIN_WISE = "WISE"
TURIN_NAME_SOURCE = "IBIS4CAT-IGR"
TURIN_RAh = "RAh"
TURIN_RAm = "RAm"
TURIN_RAs = "RAs"
TURIN_DEd = "DEd"
TURIN_DEm = "DEm"
TURIN_DEs = "DEs"
TURIN_z = "z"
TURIN_HEADER = [
    TURIN_ID,
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
TURIN_NAMES_HEADER = [
    TURIN_ID,
    "SYCAT",
    "SUMSS",
    "NVSS",
    TURIN_WISE,
    "Pan-STARRS",
    "ROSAT",
    "3PBC",
    TURIN_SWIFT,
    TURIN_NAME_SOURCE,
]
TURIN_NAMES_COLSPECS = [
    (0, 2),
    (4, 17),
    (19, 38),
    (40, 58),
    (60, 79),
    (80, 105),
    (107, 122),
    (124, 140),
    (142, 161),
    (162, 182),
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
MJD_04102023 = 60132
MJD_GOLDBRONZE_START = 58635
DAYS_IN_YEAR = 365
SECONDS_IN_YEAR = DAYS_IN_YEAR * 24 * 60 * 60
FLUX_NU = 1e36  # s-1 GeV-1 at 100 GeV
E0 = 100  # GeV
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
RATIO_90_TO_SIGMA = 2.146
RATIO_68_TO_SIGMA = 1.515
RATIO_50_TO_SIGMA = 1.177
TOTAL_SCRAMBLINGS_SPLINEMPE_TURIN = 200000
TOTAL_SCRAMBLINGS_SPLINEMPE_MILLIQUAS = 50000
TOTAL_SCRAMBLINGS_MILLIPEDE_TURIN = 1000
TOTAL_SCRAMBLINGS_MILLIPEDE_MILLIQUAS = 100
ROUND_ANGLE = 360  # deg
SPLINEMPE_ANG_DIST_FAST_SELECTION = 4  # deg
MILLIPEDE_ANG_DIST_FAST_SELECTION = 5  # deg
SPLINEMPE_SEARCH_RADIUS = 3  # deg
MILLIPEDE_SEARCH_RADIUS = 4  # deg
TEST_STATISTIC_EMPTY_SCRAMBLE = -1000

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

if catalog == ALLOWED_CATALOGS[TURIN_INDEX]:
    names_filename = TURIN_NAMES_FILENAME
    names_path = data_path / names_filename

    print(f"Checking if '{names_filename}' is in '{data_path}'...")

    if os.path.isfile(names_path):
        print(f"'{names_filename}' in '{data_path}', no need to download")
    else:
        url_names = TURIN_NAMES_URL
        print(f"{names_filename} not found, download from {url_names}...")
        r_names = requests.get(url_names, allow_redirects=True)
        open(names_path, "wb").write(r_names.content)

print(f"Loading the {catalog} catalog...")

if catalog == ALLOWED_CATALOGS[TURIN_INDEX]:
    dataframe = pd.read_fwf(
        data_path / TURIN_FILENAME, header=TURIN_HEADER_PRESENT, names=TURIN_HEADER
    )
    dataframe_names = pd.read_fwf(
        names_path,
        colspecs=TURIN_NAMES_COLSPECS,
        header=TURIN_HEADER_PRESENT,
        names=TURIN_NAMES_HEADER,
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
    names_catalog = dataframe_names[TURIN_NAME_SOURCE].to_numpy()
    names_catalog[pd.isna(names_catalog)] = dataframe_names[TURIN_SWIFT][
        pd.isna(names_catalog)
    ]
    names_catalog[pd.isna(names_catalog)] = dataframe_names[TURIN_WISE][
        pd.isna(names_catalog)
    ]
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
    names_catalog = dataframe_nearby_no0[MILLIQUAS_NAME].to_numpy()

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
days = MJD_04102023 - MJD_GOLDBRONZE_START
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


def angular_dist_score(az_true, zen_true, az_pred, zen_pred):
    """
    calculate the MAE of the angular distance between two directions.
    The two vectors are first converted to cartesian unit vectors,
    and then their scalar product is computed, which is equal to
    the cosine of the angle between the two vectors. The inverse
    cosine (arccos) thereof is then the angle between the two input vectors

    Parameters:
    -----------

    az_true : float (or array thereof)
        true azimuth value(s) in radian
    zen_true : float (or array thereof)
        true zenith value(s) in radian
    az_pred : float (or array thereof)
        predicted azimuth value(s) in radian
    zen_pred : float (or array thereof)
        predicted zenith value(s) in radian

    Returns:
    --------

    dist : float
        mean over the angular distance(s) in radian
    """
    if not (
        np.all(np.isfinite(az_true))
        and np.all(np.isfinite(zen_true))
        and np.all(np.isfinite(az_pred))
        and np.all(np.isfinite(zen_pred))
    ):
        raise ValueError("All arguments must be finite")
    # pre-compute all sine and cosine values
    sa1 = np.sin(az_true)
    ca1 = np.cos(az_true)
    sz1 = np.sin(zen_true)
    cz1 = np.cos(zen_true)
    sa2 = np.sin(az_pred)
    ca2 = np.cos(az_pred)
    sz2 = np.sin(zen_pred)
    cz2 = np.cos(zen_pred)
    # scalar product of the two cartesian vectors (x = sz*ca, y = sz*sa, z = cz)
    scalar_prod = sz1 * sz2 * (ca1 * ca2 + sa1 * sa2) + (cz1 * cz2)
    # scalar product of two unit vectors is always between -1 and 1, this is against nummerical instability
    # that might otherwise occure from the finite precision of the sine and cosine functions
    scalar_prod = np.clip(scalar_prod, -1, 1)
    # convert back to an angle (in radian)
    return np.average(np.abs(np.arccos(scalar_prod)))


def dir_contribute(ra1, dec1, ra2, dec2, raS, decS, sigma1, sigma2):  # in radiants
    """
    Contribute to the test statistic related to
    how near are the alerts to the source. All
    directions must be given in radiants
    """
    phi1 = angular_dist_score(ra1, dec1 + np.pi / 2.0, raS, decS + np.pi / 2.0)
    phi2 = angular_dist_score(ra2, dec2 + np.pi / 2.0, raS, decS + np.pi / 2.0)
    cont = -0.5 * ((phi1 / sigma1) ** 2 + (phi2 / sigma2) ** 2)
    return cont


def noise_contribute(dec1, dec2):
    """
    Contribute to the test statistic related to
    the null hypothesis
    """
    return -np.log(np.cos(dec1) * np.cos(dec2))


def test_statistic(ra1, dec1, sigma1, ra2, dec2, sigma2, raS, decS, z):
    """
    Test statistic related to two alerts and a source
    """
    c1 = flux_contribute(z, decS)
    c2 = unc_contribute(sigma1, sigma2)
    c3 = dir_contribute(ra1, dec1, ra2, dec2, raS, decS, sigma1, sigma2)
    c4 = noise_contribute(dec1, dec2)
    ts = c1 + c2 + c3 + c4
    return ts


print(f"Retrieving the alerts reconstructed with {reco}...")

if reco == ALLOWED_RECONSTRUCTIONS[SPLINEMPE_INDEX]:
    RAs = np.array([])
    DECs = np.array([])
    sigmas = np.array([])
    NAMEs = np.array([])
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
                    if len(NAMEs) > 0:
                        if name in NAMEs:
                            NAMEs[np.where(NAMEs == name)] = f"IC{year}{month}{day}B"
                    if (
                        name in rev1_names or name in SPLINEMPE_EXCEPTIONS
                    ) and not name in SPLINEMPE_BACKGROUND:
                        NAMEs = np.append(NAMEs, name)
                        has_rev1 = True
                    else:
                        has_rev1 = False
            if rev0 and is_data and has_rev1:
                if notice_line_index == 7:
                    ra = float(line.split(">")[1].split("<")[0])
                    RAs = np.append(RAs, ra)
                if notice_line_index == 8:
                    dec = float(line.split(">")[1].split("<")[0])
                    DECs = np.append(DECs, dec)
                if notice_line_index == 10:
                    err_50 = float(line.split(">")[1].split("<")[0]) / 60
                    sigma = np.deg2rad(err_50) / RATIO_50_TO_SIGMA
                    sigmas = np.append(sigmas, sigma)
elif reco == ALLOWED_RECONSTRUCTIONS[MILLIPEDE_INDEX]:
    alerts_df = pd.read_csv(data_path / MILLIPEDE_FILENAME)
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

    sigmas = np.array([])
    for i in range(len(alerts_df)):
        area = millipede_area(i)
        sigma = np.sqrt(area / np.pi) / RATIO_90_TO_SIGMA
        sigmas = np.append(sigmas, sigma)

print("Estimating background...")

if reco == ALLOWED_RECONSTRUCTIONS[SPLINEMPE_INDEX]:
    ang_dist_fast_selection = SPLINEMPE_ANG_DIST_FAST_SELECTION
    search_radius = SPLINEMPE_SEARCH_RADIUS
    if catalog == ALLOWED_CATALOGS[TURIN_INDEX]:
        total_scramblings = TOTAL_SCRAMBLINGS_SPLINEMPE_TURIN
    elif catalog == ALLOWED_CATALOGS[MILLIQUAS_INDEX]:
        total_scramblings = TOTAL_SCRAMBLINGS_SPLINEMPE_MILLIQUAS
elif reco == ALLOWED_RECONSTRUCTIONS[MILLIPEDE_INDEX]:
    ang_dist_fast_selection = MILLIPEDE_ANG_DIST_FAST_SELECTION
    search_radius = MILLIPEDE_SEARCH_RADIUS
    if catalog == ALLOWED_CATALOGS[TURIN_INDEX]:
        total_scramblings = TOTAL_SCRAMBLINGS_MILLIPEDE_TURIN
    elif catalog == ALLOWED_CATALOGS[MILLIQUAS_INDEX]:
        total_scramblings = TOTAL_SCRAMBLINGS_MILLIPEDE_MILLIQUAS
test_statistic_per_scramble = np.array([])
names_alerts_per_scramble = np.array([])
names_source_per_scramble = np.array([])
t0 = time.time()
for scrambling_number in range(total_scramblings):
    if (scrambling_number+1)%100 == 0:
        sys.stdout.write(
            "\r"
            + f"Scrumble nr {scrambling_number + 1:6}"
            + " of "
            + str(total_scramblings)
            + f". Taken {round(time.time() - t0, 1):6}"
            + f" seconds so far. Still {round((total_scramblings - (scrambling_number + 1)) * (time.time() - t0) / (scrambling_number + 1), 1):6}"
            + " seconds remaining."
        )
    rng = np.random.default_rng(seed=scrambling_number)
    random_ras = rng.uniform(0.0, ROUND_ANGLE, size=len(RAs))
    test_statistic_per_doublet = np.array([])
    names_alerts_per_doublet = np.array([])
    names_source_per_doublet = np.array([])
    for first_alert_index in range(len(RAs)):
        for second_alert_index in range(first_alert_index, len(RAs)):
            if first_alert_index == second_alert_index:
                continue
            ra1 = random_ras[first_alert_index]
            ra2 = random_ras[second_alert_index]
            dec1 = DECs[first_alert_index]
            dec2 = DECs[second_alert_index]
            higher_ra = max(ra1, ra2)
            smaller_ra = min(ra1, ra2)
            # First fast selection
            ra_distance = min(higher_ra - smaller_ra, smaller_ra - higher_ra + 360)
            dec_distance = np.abs(dec1 - dec2)
            if (
                ra_distance > ang_dist_fast_selection
                or dec_distance > ang_dist_fast_selection
            ):
                continue
            ra1_rad = np.deg2rad(ra1)
            ra2_rad = np.deg2rad(ra2)
            dec1_rad = np.deg2rad(dec1)
            dec2_rad = np.deg2rad(dec2)
            sigma1 = sigmas[first_alert_index]
            sigma2 = sigmas[second_alert_index]
            # Consider the sources nearest to the alert with best angular resolution
            if sigma1 <= sigma2:
                search_index = first_alert_index
            else:
                search_index = second_alert_index
            mask_ra = np.logical_and(
                RAs_catalog < random_ras[search_index] + search_radius,
                RAs_catalog > random_ras[search_index] - search_radius,
            )
            mask_dec = np.logical_and(
                DECs_catalog < DECs[search_index] + search_radius,
                DECs_catalog > DECs[search_index] - search_radius,
            )
            mask_sources = np.logical_and(mask_ra, mask_dec)
            test_statistic_per_source = np.array([])
            RAs_sources_nearby = RAs_catalog[mask_sources]
            DECs_sources_nearby = DECs_catalog[mask_sources]
            redshifts_sources_nearby = redshifts_catalog[mask_sources]
            names_sources_nearby = names_catalog[mask_sources]
            for source_index in range(len(names_sources_nearby)):
                ra_rad_source = np.deg2rad(RAs_sources_nearby[source_index])
                de_rad_source = np.deg2rad(DECs_sources_nearby[source_index])
                redshift_source = redshifts_sources_nearby[source_index]
                test_statistic_source = test_statistic(
                    ra1_rad,
                    dec1_rad,
                    sigma1,
                    ra2_rad,
                    dec2_rad,
                    sigma2,
                    ra_rad_source,
                    de_rad_source,
                    redshift_source,
                )
                test_statistic_per_source = np.append(
                    test_statistic_per_source, test_statistic_source
                )
            if len(test_statistic_per_source) == 0:
                continue
            index_best_ts = np.argmax(test_statistic_per_source)
            test_statistic_doublet = test_statistic_per_source[index_best_ts]
            name_best_source = names_sources_nearby[index_best_ts]
            test_statistic_per_doublet = np.append(
                test_statistic_per_doublet, test_statistic_doublet
            )
            names_alerts_per_doublet = np.append(
                names_alerts_per_doublet,
                f"{NAMEs[first_alert_index]}, {NAMEs[second_alert_index]}",
            )
            names_source_per_doublet = np.append(
                names_source_per_doublet, name_best_source
            )
    if len(test_statistic_per_doublet) == 0:
        test_statistic_scramble = TEST_STATISTIC_EMPTY_SCRAMBLE
        names_alerts_scramble = None
        name_source_scramble = None
    else:
        index_best_doublet = np.argmax(test_statistic_per_doublet)
        test_statistic_scramble = test_statistic_per_doublet[index_best_doublet]
        names_alerts_scramble = names_alerts_per_doublet[index_best_doublet]
        name_source_scramble = names_source_per_doublet[index_best_doublet]
    test_statistic_per_scramble = np.append(
        test_statistic_per_scramble, test_statistic_scramble
    )
    names_alerts_per_scramble = np.append(
        names_alerts_per_scramble, names_alerts_scramble
    )
    names_source_per_scramble = np.append(
        names_source_per_scramble, name_source_scramble
    )

print(f"\nEstimate ts value for {reco} with {catalog}...")

test_statistic_per_doublet = np.array([])
names_alerts_per_doublet = np.array([])
names_source_per_doublet = np.array([])
for first_alert_index in range(len(RAs)):
    for second_alert_index in range(first_alert_index, len(RAs)):
        if first_alert_index == second_alert_index:
            continue
        ra1 = RAs[first_alert_index]
        ra2 = RAs[second_alert_index]
        dec1 = DECs[first_alert_index]
        dec2 = DECs[second_alert_index]
        higher_ra = max(ra1, ra2)
        smaller_ra = min(ra1, ra2)
        # First fast selection
        ra_distance = min(higher_ra - smaller_ra, smaller_ra - higher_ra + 360)
        dec_distance = np.abs(dec1 - dec2)
        if (
            ra_distance > ang_dist_fast_selection
            or dec_distance > ang_dist_fast_selection
        ):
            continue
        ra1_rad = np.deg2rad(ra1)
        ra2_rad = np.deg2rad(ra2)
        dec1_rad = np.deg2rad(dec1)
        dec2_rad = np.deg2rad(dec2)
        sigma1 = sigmas[first_alert_index]
        sigma2 = sigmas[second_alert_index]
        # Consider the sources nearest to the alert with best angular resolution
        if sigma1 <= sigma2:
            search_index = first_alert_index
        else:
            search_index = second_alert_index
        mask_ra = np.logical_and(
            RAs_catalog < RAs[search_index] + search_radius,
            RAs_catalog > RAs[search_index] - search_radius,
        )
        mask_dec = np.logical_and(
            DECs_catalog < DECs[search_index] + search_radius,
            DECs_catalog > DECs[search_index] - search_radius,
        )
        mask_sources = np.logical_and(mask_ra, mask_dec)
        test_statistic_per_source = np.array([])
        RAs_sources_nearby = RAs_catalog[mask_sources]
        DECs_sources_nearby = DECs_catalog[mask_sources]
        redshifts_sources_nearby = redshifts_catalog[mask_sources]
        names_sources_nearby = names_catalog[mask_sources]
        for source_index in range(len(names_sources_nearby)):
            ra_rad_source = np.deg2rad(RAs_sources_nearby[source_index])
            de_rad_source = np.deg2rad(DECs_sources_nearby[source_index])
            redshift_source = redshifts_sources_nearby[source_index]
            test_statistic_source = test_statistic(
                ra1_rad,
                dec1_rad,
                sigma1,
                ra2_rad,
                dec2_rad,
                sigma2,
                ra_rad_source,
                de_rad_source,
                redshift_source,
            )
            test_statistic_per_source = np.append(
                test_statistic_per_source, test_statistic_source
            )
        if len(test_statistic_per_source) == 0:
            continue
        index_best_ts = np.argmax(test_statistic_per_source)
        test_statistic_doublet = test_statistic_per_source[index_best_ts]
        name_best_source = names_sources_nearby[index_best_ts]
        test_statistic_per_doublet = np.append(
            test_statistic_per_doublet, test_statistic_doublet
        )
        names_alerts_per_doublet = np.append(
            names_alerts_per_doublet,
            f"{NAMEs[first_alert_index]}, {NAMEs[second_alert_index]}",
        )
        names_source_per_doublet = np.append(names_source_per_doublet, name_best_source)
best_test_statistic_index = np.argmax(test_statistic_per_doublet)
best_test_statistic = test_statistic_per_doublet[best_test_statistic_index]
best_source = names_source_per_doublet[best_test_statistic_index]
best_alerts = names_alerts_per_doublet[best_test_statistic_index]

print(
    f"\nTest statistic for {reco} with {catalog}: {best_test_statistic}\nSource: {best_source}. Doublet: {best_alerts}"
)

p_value = len(
    test_statistic_per_scramble[test_statistic_per_scramble > best_test_statistic]
) / len(test_statistic_per_scramble)

print(
    f"\n\n*********************\n\n"
    f"p-value: {p_value} = {p_value*100}%"
    f"\n\n*********************\n\n"
)
