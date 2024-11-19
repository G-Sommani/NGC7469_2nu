from abc import ABC, abstractmethod
import importlib
import numpy as np
from pathlib import Path
from typing import List
from scipy.stats import multinomial # type: ignore


class Catalog(ABC):

    catalog_name: str = __name__
    xray: bool = False
    noweight: bool = False
    zipname_data: str | None = None
    filename_names: str | None = None
    url_names: str | None = None
    filename_xray: str | None = None
    url_xray: str | None = None

    ras_catalog: np.ndarray = np.array([])
    decs_catalog: np.ndarray = np.array([])
    redshifts_catalog: np.ndarray = np.array([])
    names_catalog: np.ndarray = np.array([])
    xray_catalog: np.ndarray = np.array([])

    filename_data: str
    url_data: str
    total_scrambling_possibilities: List[int]

    def __init__(
            self, xray: bool = False, noweight: bool = False,
    ) -> None:
        self.xray = xray
        self.noweight = noweight
        pass

    @abstractmethod
    def load_catalog(self, data_path: Path) -> None:
        pass


def initiate_catalog(
        catalog_name: str, xray: bool = False, noweight: bool = False
    ) -> Catalog:
    module = importlib.import_module(f"{__name__}.{catalog_name.lower()}")
    return module.CATALOG_CLASS(xray=xray, noweight=noweight)
