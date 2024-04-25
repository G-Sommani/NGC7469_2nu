import os
from pathlib import Path
from typing import Tuple
import config as cfg
import zipfile
import requests


def define_catalog(catalog: str) -> Tuple[str, str]:

    filename_data = None
    url_data = None
    filename_names = None
    url_names = None

    if catalog == cfg.ALLOWED_CATALOGS[cfg.TURIN_INDEX]:
        filename_data = cfg.TURIN_FILENAME
        url_data = cfg.TURIN_URL
        filename_names = cfg.TURIN_NAMES_FILENAME
        url_names = cfg.TURIN_NAMES_URL
    elif catalog == cfg.ALLOWED_CATALOGS[cfg.MILLIQUAS_INDEX]:
        filename_data = cfg.MILLIQUAS_FILENAME
        url_data = cfg.MILLIQUAS_URL

    return filename_data, url_data, filename_names, url_names


class Loader:
    def __init__(self):

        print("Definition of paths...")

        # Definition of paths
        self.cwd = Path(os.getcwd())
        self.data_path = self.cwd / "../data"
        self.data_results_path = self.cwd / "../data_results"
        self.figures_path = self.cwd / "../figures"

    def download_catalog(self, catalog: str) -> None:

        filename_data, url_data, filename_names, url_names = define_catalog(catalog)

        print(f"Checking if '{filename_data}' is in '{self.data_path}'...")

        if os.path.isfile(self.data_path / filename_data):
            print(f"'{filename_data}' in '{self.data_path}', no need to download")

        else:
            print(f"{filename_data} not found, download {catalog} catalog...")

            r = requests.get(url_data, allow_redirects=True)
            if catalog == cfg.ALLOWED_CATALOGS[cfg.TURIN_INDEX]:
                catalog_path = self.data_path / cfg.TURIN_FILENAME
                open(catalog_path, "wb").write(r.content)

                print(f"Checking if '{filename_names}' is in '{self.data_path}'...")

                path_names = self.data_path / filename_names

                if os.path.isfile(path_names):
                    print(
                        f"'{filename_names}' in '{self.data_path}', no need to download"
                    )

                else:
                    print(f"{filename_names} not found, download from {url_names}...")
                    r_names = requests.get(url_names, allow_redirects=True)
                    open(path_names, "wb").write(r_names.content)

            elif catalog == cfg.ALLOWED_CATALOGS[cfg.MILLIQUAS_INDEX]:
                zip_path = self.data_path / cfg.MILLIQUAS_ZIP
                open(zip_path, "wb").write(r.content)

                print(f"Unzipping {catalog} catalog...")

                with zipfile.ZipFile(zip_path, "r") as zip_ref:
                    zip_ref.extractall(self.data_path)

                print(f"Removing {cfg.MILLIQUAS_ZIP}...")

                os.remove(zip_path)

        return
