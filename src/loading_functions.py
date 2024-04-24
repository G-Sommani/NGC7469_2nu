import os
from pathlib import Path
from typing import Tuple
import config as cfg
import zipfile

def define_catalog(catalog: str) -> Tuple[str, str]:

    filename = None
    url = None
    if catalog == cfg.ALLOWED_CATALOGS[cfg.TURIN_INDEX]:
        filename = cfg.TURIN_FILENAME
        url = cfg.TURIN_URL
    elif catalog == cfg.ALLOWED_CATALOGS[cfg.MILLIQUAS_INDEX]:
        filename = cfg.MILLIQUAS_FILENAME
        url = cfg.MILLIQUAS_URL

    return filename, url

class Loader:

    def __init__(self):

        print("Definition of paths...")

        # Definition of paths
        self.cwd = Path(os.getcwd())
        self.data_path = self.cwd / "../data"
        self.data_results_path = self.cwd / "../data_results"
        self.figures_path = self.cwd / "../figures"

    def download_catalog(self, catalog: str) -> None:

        filename, url = define_catalog(catalog)

        print(f"Checking if '{filename}' is in '{self.data_path}'...")

        if os.path.isfile(self.data_path / filename):
            print(f"'{filename}' in '{self.data_path}', no need to download")

        else:
            print(f"{filename} not found, download {catalog} catalog...")

            r = requests.get(url, allow_redirects=True)
            if catalog == cfg.ALLOWED_CATALOGS[cfg.TURIN_INDEX]:
                catalog_path = self.data_path / cfg.TURIN_FILENAME
                open(catalog_path, "wb").write(r.content)
            elif catalog == cfg.ALLOWED_CATALOGS[cfg.MILLIQUAS_INDEX]:
                zip_path = self.data_path / cfg.MILLIQUAS_ZIP
                open(zip_path, "wb").write(r.content)

                print(f"Unzipping {catalog} catalog...")

                with zipfile.ZipFile(zip_path, "r") as zip_ref:
                    zip_ref.extractall(self.data_path)

                print(f"Removing {cfg.MILLIQUAS_ZIP}...")

                os.remove(zip_path)

        return