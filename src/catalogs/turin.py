from typing import Final, List
import pandas as pd  # type: ignore
from pathlib import Path
import numpy as np

from . import Catalog
from . import config_catalogs as cfg


class Turin(Catalog):
    def __init__(self, xray: bool = False, noweight: bool = False) -> None:
        super().__init__(xray=xray)
        self.catalog_name = cfg.TURIN
        self.filename_data = cfg.TURIN_FILENAME
        self.url_data = cfg.TURIN_URL
        self.filename_names = cfg.TURIN_NAMES_FILENAME
        self.url_names = cfg.TURIN_NAMES_URL
        if self.xray:
            self.filename_xray = cfg.BASS_XRAY_FILENAME
            self.url_xray = cfg.BASS_XRAY_URL
            self.total_scrambling_possibilities = [
                cfg.TOTAL_SCRAMBLINGS_SPLINEMPE_TURIN_XRAY,
                cfg.TOTAL_SCRAMBLINGS_MILLIPEDE_TURIN_XRAY,
            ]
        else:
            self.total_scrambling_possibilities = [
                cfg.TOTAL_SCRAMBLINGS_SPLINEMPE_TURIN,
                cfg.TOTAL_SCRAMBLINGS_SPLINEMPE_TURIN,
            ]

    def load_catalog(self, data_path: Path) -> None:

        dataframe = pd.read_fwf(
            data_path / cfg.TURIN_FILENAME,
            header=cfg.TURIN_HEADER_PRESENT,
            names=cfg.TURIN_HEADER,
        )
        dataframe_names = pd.read_fwf(
            data_path / cfg.TURIN_NAMES_FILENAME,
            colspecs=cfg.TURIN_NAMES_COLSPECS,
            header=cfg.TURIN_HEADER_PRESENT,
            names=cfg.TURIN_NAMES_HEADER,
        )
        if self.xray:
            dataframe_flux = pd.read_csv(
                data_path / cfg.BASS_XRAY_FILENAME,
                names=cfg.BASS_XRAY_HEADER,
                skiprows=cfg.BASS_SKIPROWS,
                sep=cfg.BASS_SEP,
            )

            # Drop sources in no BAT-Swift catalog.
            turin_swift_names_array = dataframe_names[cfg.TURIN_SWIFT].to_numpy()
            turin_swift_names_list = list(turin_swift_names_array)
            (idxs_sources_no_intr,) = np.where(pd.isnull(turin_swift_names_array))
            sycat_names = dataframe_names[cfg.TURIN_SYCAT].to_numpy()
            print(
                "Dropping the following sources (SYCAT IDs) because they are not in any BAT Swift catalog:"
            )
            print(sycat_names[idxs_sources_no_intr])
            dataframe = dataframe.drop(idxs_sources_no_intr)
            dataframe_names = dataframe_names.drop(idxs_sources_no_intr)

            # Drop sources that are in Swift BAT 105-Month, but not in Swift BAT 70-Month
            swift_names = dataframe_names[cfg.TURIN_SWIFT].to_numpy()
            flux_source_names = list(dataframe_flux[cfg.BASS_SWIFT].to_numpy())
            notfound_sources = list()
            notfound_idxs = list()
            for i, name in enumerate(swift_names):
                namecode = name[:5] + name[6:]
                if namecode not in flux_source_names:
                    notfound_sources.append(name)
                    idx = turin_swift_names_list.index(name)
                    notfound_idxs.append(idx)
            print(
                "Dropping the following sources (Swift BAT IDs) because they are not in the 70 Month BAT Swift catalog:"
            )
            print(turin_swift_names_array[notfound_idxs])
            dataframe = dataframe.drop(notfound_idxs)
            dataframe_names = dataframe_names.drop(notfound_idxs)

            # Define array with intrinsic x-ray fluxes
            swift_bat_names = dataframe_names[cfg.TURIN_SWIFT].to_numpy()
            all_xray_fluxes = dataframe_flux[cfg.XRAY_FLUX_NAME].to_numpy()
            for name in swift_bat_names:
                namecode = name[:5] + name[6:]
                name_idx = flux_source_names.index(namecode)
                intr_xray_flux = all_xray_fluxes[name_idx]
                self.xray_catalog = np.append(
                    self.xray_catalog, intr_xray_flux * cfg.FLUX_FACTOR
                )

        def hms_to_deg(h, m, s):
            return 15.0 * (h + (m + s / 60.0) / 60.0)

        def dms_to_deg(d, m, s):
            return d + (m + s / 60.0) / 60.0

        self.ras_catalog = hms_to_deg(
            dataframe[cfg.TURIN_RAh],
            dataframe[cfg.TURIN_RAm],
            dataframe[cfg.TURIN_RAs],
        ).to_numpy()
        self.decs_catalog = dms_to_deg(
            dataframe[cfg.TURIN_DEd],
            dataframe[cfg.TURIN_DEm],
            dataframe[cfg.TURIN_DEs],
        ).to_numpy()
        self.redshifts_catalog = dataframe[cfg.TURIN_z].to_numpy()
        self.names_catalog = dataframe_names[cfg.TURIN_NAME_SOURCE].to_numpy()
        self.names_catalog[pd.isna(self.names_catalog)] = dataframe_names[
            cfg.TURIN_SWIFT
        ][pd.isna(self.names_catalog)]
        self.names_catalog[pd.isna(self.names_catalog)] = dataframe_names[
            cfg.TURIN_WISE
        ][pd.isna(self.names_catalog)]


CATALOG_CLASS: Final[type[Catalog]] = Turin
