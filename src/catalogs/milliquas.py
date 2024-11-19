from typing import Final, List
import pandas as pd  # type: ignore
import numpy as np
from pathlib import Path

from . import Catalog
from . import config_catalogs as cfg


def agn_mask(dataframe: pd.DataFrame) -> np.ndarray:

    mask01 = np.ma.mask_or(
        dataframe[cfg.MILLIQUAS_TYPE] == cfg.MILLIQUAS_AGN_CATEGORIES[0],
        dataframe[cfg.MILLIQUAS_TYPE] == cfg.MILLIQUAS_AGN_CATEGORIES[1],
    )
    mask23 = np.ma.mask_or(
        dataframe[cfg.MILLIQUAS_TYPE] == cfg.MILLIQUAS_AGN_CATEGORIES[2],
        dataframe[cfg.MILLIQUAS_TYPE] == cfg.MILLIQUAS_AGN_CATEGORIES[3],
    )
    mask45 = np.ma.mask_or(
        dataframe[cfg.MILLIQUAS_TYPE] == cfg.MILLIQUAS_AGN_CATEGORIES[4],
        dataframe[cfg.MILLIQUAS_TYPE] == cfg.MILLIQUAS_AGN_CATEGORIES[5],
    )
    mask67 = np.ma.mask_or(
        dataframe[cfg.MILLIQUAS_TYPE] == cfg.MILLIQUAS_AGN_CATEGORIES[6],
        dataframe[cfg.MILLIQUAS_TYPE] == cfg.MILLIQUAS_AGN_CATEGORIES[7],
    )
    mask89 = np.ma.mask_or(
        dataframe[cfg.MILLIQUAS_TYPE] == cfg.MILLIQUAS_AGN_CATEGORIES[8],
        dataframe[cfg.MILLIQUAS_TYPE] == cfg.MILLIQUAS_AGN_CATEGORIES[9],
    )
    mask0123 = np.ma.mask_or(mask01, mask23)
    mask4567 = np.ma.mask_or(mask45, mask67)
    mask01234567 = np.ma.mask_or(mask0123, mask4567)
    mask_agn = np.ma.mask_or(mask01234567, mask89)

    return mask_agn


class Milliquas(Catalog):
    def __init__(self, xray: bool = False, noweight: bool = False) -> None:
        if xray:
            raise ValueError(
                f"Possible to use the x-ray fluxes as weighting only with the Turin catalog."
            )
        super().__init__()
        self.catalog_name = cfg.MILLIQUAS
        self.filename_data = cfg.MILLIQUAS_FILENAME
        self.url_data = cfg.MILLIQUAS_URL
        self.zipname_data = cfg.MILLIQUAS_ZIP
        self.noweight = noweight
        if self.noweight:
            total_scrambling_splinempe = cfg.TOTAL_SCRAMBLINGS_SPLINEMPE_MILLIQUAS_NOWEIGHT
        else:
            total_scrambling_splinempe = cfg.TOTAL_SCRAMBLINGS_SPLINEMPE_MILLIQUAS_WEIGHT
        self.total_scrambling_possibilities = [
            total_scrambling_splinempe,
            cfg.TOTAL_SCRAMBLINGS_MILLIPEDE_MILLIQUAS,
        ]

    def mask_sources(self, dataframe: pd.DataFrame) -> pd.DataFrame:

        # Consider only the agn in the milliquas catalog
        mask_agn = agn_mask(dataframe)
        dataframe_agn = dataframe[mask_agn]
        print(f"noweight is {self.noweight}")
        if not self.noweight:
            print(f"Selecting only AGN within a redshift of {cfg.MILLIQUAS_Z_CUT}..")
            # Sources that are too far away are not significant for the test and slow down the program
            dataframe_agn = dataframe_agn[
                dataframe_agn[cfg.MILLIQUAS_Z] < cfg.MILLIQUAS_Z_CUT
            ]

        # There are two blazars which are reported with a redshift of zero. This is unphyisical and
        # these two sources are therefore removed.
        dataframe_nearby_no0 = dataframe_agn[dataframe_agn[cfg.MILLIQUAS_Z] != 0]

        return dataframe_nearby_no0

    def load_catalog(self, data_path: Path) -> None:

        dataframe = pd.read_fwf(
            data_path / cfg.MILLIQUAS_FILENAME,
            colspecs=cfg.MILLIQUAS_COLSPECS,
            names=cfg.MILLIQUAS_HEADER,
        )

        masked_dataframe = self.mask_sources(dataframe)

        self.ras_catalog = masked_dataframe[cfg.MILLIQUAS_RA].to_numpy()
        self.decs_catalog = masked_dataframe[cfg.MILLIQUAS_DEC].to_numpy()
        self.redshifts_catalog = masked_dataframe[cfg.MILLIQUAS_Z].to_numpy()
        self.names_catalog = masked_dataframe[cfg.MILLIQUAS_NAME].to_numpy()


CATALOG_CLASS: Final[type[Catalog]] = Milliquas
