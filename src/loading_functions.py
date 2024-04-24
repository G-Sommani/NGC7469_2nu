import os
from pathlib import Path
from typing import List, Tuple
import config as cfg

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