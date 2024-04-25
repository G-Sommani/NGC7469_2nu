from typing import Final

from . import Catalog
from . import config_catalogs as cfg


class Turin(Catalog):
    def __init__(self, xray: bool = False) -> None:
        super().__init__(xray=xray)
        self.filename_data = cfg.TURIN_FILENAME
        self.url_data = cfg.TURIN_URL
        self.filename_names = cfg.TURIN_NAMES_FILENAME
        self.url_names = cfg.TURIN_NAMES_URL
        if self.xray:
            self.filename_xray = cfg.BASS_XRAY_FILENAME
            self.url_xray = cfg.BASS_XRAY_URL


CATALOG_CLASS: Final[type[Catalog]] = Turin
