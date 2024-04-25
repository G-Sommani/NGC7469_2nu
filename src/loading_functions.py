import os
from pathlib import Path
import zipfile
import requests
import catalogs
import numpy as np
import config as cfg


class Loader:
    def __init__(self):

        print("Definition of paths...")

        # Definition of paths
        self.cwd = Path(os.getcwd())
        self.data_path = self.cwd / "../data"
        self.data_results_path = self.cwd / "../data_results"
        self.figures_path = self.cwd / "../figures"

    def download_file(
        self, filename: str | None, url: str | None, zipname: str | None = None
    ) -> None:
        print(f"{filename} not found, download from {url}...")
        r = requests.get(str(url), allow_redirects=True)
        if zipname is not None:
            zip_path = self.data_path / zipname
            open(zip_path, "wb").write(r.content)

            print(f"Unzipping {zipname} catalog...")

            with zipfile.ZipFile(zip_path, "r") as zip_ref:
                zip_ref.extractall(self.data_path)

            print(f"Removing {zipname}...")

            os.remove(zip_path)
        else:
            open(self.data_path / filename, "wb").write(r.content)

    def check_presence(
        self, filename: str | None, url: str | None, zipname: str | None = None
    ) -> None:

        if filename == None:
            return

        print(f"Checking if '{filename}' is in '{self.data_path}'...")

        if os.path.isfile(self.data_path / filename):
            print(f"'{filename}' in '{self.data_path}', no need to download")

        else:
            self.download_file(filename, url, zipname=zipname)

    def download_catalog(self, catalog: catalogs.Catalog) -> None:
        self.check_presence(
            catalog.filename_data, catalog.url_data, zipname=catalog.zipname_data
        )
        self.check_presence(catalog.filename_names, catalog.url_names)
        self.check_presence(catalog.filename_xray, catalog.url_xray)

    def load_catalog(self, catalog: catalogs.Catalog) -> None:

        self.download_catalog(catalog)

        print(f"Loading the {catalog.catalog_name} catalog...")

        catalog.load_catalog(self.data_path)

    def load_effective_area(self) -> np.ndarray:
        return np.genfromtxt(self.data_path / cfg.EFFECTIVE_AREA_FILENAME)
