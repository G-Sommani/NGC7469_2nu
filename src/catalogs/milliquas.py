from typing import Final

from . import Catalog
from . import config_catalogs as cfg


class Milliquas(Catalog):
    def __init__(self, xray: bool = False) -> None:
        if xray:
            raise ValueError(
                f"Possible to use the x-ray fluxes as weighting only with the Turin catalog."
            )
        super().__init__()
        self.filename_data = cfg.MILLIQUAS_FILENAME
        self.url_data = cfg.MILLIQUAS_URL
        self.zipname_data = cfg.MILLIQUAS_ZIP


CATALOG_CLASS: Final[type[Catalog]] = Milliquas
