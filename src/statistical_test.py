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
    print(len(dataframe_nearby))
    #
    dataframe_nearby_no0 = dataframe_nearby[dataframe_nearby[MILLIQUAS_Z] != 0]
    RAs_catalog = dataframe_nearby_no0[MILLIQUAS_RA].to_numpy()
    DECs_catalog = dataframe_nearby_no0[MILLIQUAS_DEC].to_numpy()
    redshifts_catalog = dataframe_nearby_no0[MILLIQUAS_Z].to_numpy()

print(redshifts_catalog)
