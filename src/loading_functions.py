import os
from pathlib import Path
from typing import List, Optional
import config as cfg
import zipfile
import requests


def define_catalog(catalog: str, flux: bool = False) -> List[Optional[str]]:

    filename_data = None
    url_data = None
    zipname_data = None
    filename_names = None
    url_names = None
    filename_xray = None
    url_xray = None

    if catalog == cfg.ALLOWED_CATALOGS[cfg.TURIN_INDEX]:
        filename_data = cfg.TURIN_FILENAME
        url_data = cfg.TURIN_URL
        filename_names = cfg.TURIN_NAMES_FILENAME
        url_names = cfg.TURIN_NAMES_URL
        if flux:
            filename_xray = cfg.BASS_XRAY_FILENAME
            url_xray = cfg.BASS_XRAY_URL
    elif catalog == cfg.ALLOWED_CATALOGS[cfg.MILLIQUAS_INDEX]:
        filename_data = cfg.MILLIQUAS_FILENAME
        url_data = cfg.MILLIQUAS_URL
        zipname_data = cfg.MILLIQUAS_ZIP

    results = [
        filename_data,
        url_data,
        zipname_data,
        filename_names,
        url_names,
        filename_xray,
        url_xray,
    ]

    return results


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

        print(f"Checking if '{filename}' is in '{self.data_path}'...")

        if os.path.isfile(self.data_path / filename):
            print(f"'{filename}' in '{self.data_path}', no need to download")

        else:
            self.download_file(filename, url, zipname=zipname)

    def download_catalog(self, catalog: str, flux: bool = False) -> None:

        (
            filename_data,
            url_data,
            zipname_data,
            filename_names,
            url_names,
            filename_xray,
            url_xray,
        ) = define_catalog(catalog, flux=flux)

        self.check_presence(filename_data, url_data, zipname=zipname_data)

        if filename_names is not None:
            self.check_presence(filename_names, url_names)

        if filename_xray is not None:
            self.check_presence(filename_xray, url_xray)
