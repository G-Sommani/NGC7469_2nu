import sys
from pathlib import Path
import os

ALLOWED_CATALOGS = ["turin", "milliquas"]
ALLOWED_RECONSTRUCTIONS = ["splinempe", "millipede"]
DEFAULT_CATALOG = ALLOWED_CATALOGS[0]
DEFAULT_RECO = ALLOWED_RECONSTRUCTIONS[0]

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
