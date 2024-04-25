import os
from pathlib import Path
from typing import List, Optional
import config as cfg
import zipfile
import requests
import numpy as np
import pandas as pd  # type: ignore
import catalogs


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

    def download_catalog(self, catalog_name: str, flux: bool = False) -> None:
        catalog = catalogs.initiate_catalog(catalog_name, xray=flux)
        self.check_presence(
            catalog.filename_data, catalog.url_data, zipname=catalog.zipname_data
        )
        self.check_presence(catalog.filename_names, catalog.url_names)
        self.check_presence(catalog.filename_xray, catalog.url_xray)

    def load_turin_catalog(self, flux: bool = False) -> List[np.ndarray]:

        ras_catalog = np.array([])
        decs_catalog = np.array([])
        redshifts_catalog = np.array([])
        names_catalog = np.array([])
        xray_catalog = np.array([])

        dataframe = pd.read_fwf(
            self.data_path / cfg.TURIN_FILENAME,
            header=cfg.TURIN_HEADER_PRESENT,
            names=cfg.TURIN_HEADER,
        )
        dataframe_names = pd.read_fwf(
            self.data_path / cfg.TURIN_NAMES_FILENAME,
            colspecs=cfg.TURIN_NAMES_COLSPECS,
            header=cfg.TURIN_HEADER_PRESENT,
            names=cfg.TURIN_NAMES_HEADER,
        )
        if flux:
            dataframe_flux = pd.read_csv(
                self.data_path / cfg.BASS_XRAY_FILENAME,
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
                xray_catalog = np.append(xray_catalog, intr_xray_flux * cfg.FLUX_FACTOR)

            def hms_to_deg(h, m, s):
                return 15.0 * (h + (m + s / 60.0) / 60.0)

            def dms_to_deg(d, m, s):
                return d + (m + s / 60.0) / 60.0

            ras_catalog = hms_to_deg(
                dataframe[cfg.TURIN_RAh],
                dataframe[cfg.TURIN_RAm],
                dataframe[cfg.TURIN_RAs],
            ).to_numpy()
            decs_catalog = dms_to_deg(
                dataframe[cfg.TURIN_DEd],
                dataframe[cfg.TURIN_DEm],
                dataframe[cfg.TURIN_DEs],
            ).to_numpy()
            redshifts_catalog = dataframe[cfg.TURIN_z].to_numpy()
            names_catalog = dataframe_names[cfg.TURIN_NAME_SOURCE].to_numpy()
            names_catalog[pd.isna(names_catalog)] = dataframe_names[cfg.TURIN_SWIFT][
                pd.isna(names_catalog)
            ]
            names_catalog[pd.isna(names_catalog)] = dataframe_names[cfg.TURIN_WISE][
                pd.isna(names_catalog)
            ]

        results = [
            ras_catalog,
            decs_catalog,
            redshifts_catalog,
            names_catalog,
            xray_catalog,
        ]

        return results

    def load_milliquas_catalog(self) -> List[np.ndarray]:

        ras_catalog = np.array([])
        decs_catalog = np.array([])
        redshifts_catalog = np.array([])
        names_catalog = np.array([])
        xray_catalog = np.array([])

        dataframe = pd.read_fwf(
            self.data_path / cfg.MILLIQUAS_FILENAME,
            colspecs=cfg.MILLIQUAS_COLSPECS,
            names=cfg.MILLIQUAS_HEADER,
        )

        print(f"Selecting only AGN within a redshift of {cfg.MILLIQUAS_Z_CUT}..")

        # Consider only the agn in the milliquas catalog
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
        dataframe_agn = dataframe[mask_agn]

        # Sources that are too far away are not significant for the test and slow down the program
        dataframe_nearby = dataframe_agn[
            dataframe_agn[cfg.MILLIQUAS_Z] < cfg.MILLIQUAS_Z_CUT
        ]

        # There are two blazars which are reported with a redshift of zero. This is unphyisical and
        # these two sources are therefore removed.
        dataframe_nearby_no0 = dataframe_nearby[dataframe_nearby[cfg.MILLIQUAS_Z] != 0]
        ras_catalog = dataframe_nearby_no0[cfg.MILLIQUAS_RA].to_numpy()
        decs_catalog = dataframe_nearby_no0[cfg.MILLIQUAS_DEC].to_numpy()
        redshifts_catalog = dataframe_nearby_no0[cfg.MILLIQUAS_Z].to_numpy()
        names_catalog = dataframe_nearby_no0[cfg.MILLIQUAS_NAME].to_numpy()

        results = [
            ras_catalog,
            decs_catalog,
            redshifts_catalog,
            names_catalog,
            xray_catalog,
        ]

        return results

    def load_catalog(self, catalog_name: str, flux: bool = False) -> List[np.ndarray]:

        self.download_catalog(catalog_name, flux=flux)

        print(f"Loading the {catalog_name} catalog...")

        if catalog_name == cfg.ALLOWED_CATALOGS[cfg.TURIN_INDEX]:
            results = self.load_turin_catalog(flux=flux)

        elif catalog_name == cfg.ALLOWED_CATALOGS[cfg.MILLIQUAS_INDEX]:
            results = self.load_milliquas_catalog()

        return results
