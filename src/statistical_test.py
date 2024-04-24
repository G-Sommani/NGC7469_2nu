import sys
import os
import requests
import zipfile
import pandas as pd
import numpy as np
import time
import argparse
import config as cfg
from loading_functions import define_paths


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--reco",
        "-r",
        type=str,
        default=cfg.ALLOWED_RECONSTRUCTIONS[cfg.SPLINEMPE_INDEX],
        help="reconstruction to use for the neutrino events",
        choices=cfg.ALLOWED_RECONSTRUCTIONS,
    )
    parser.add_argument(
        "--catalog",
        "-c",
        type=str,
        default=cfg.ALLOWED_CATALOGS[cfg.TURIN_INDEX],
        help="catalog of sources for the statistical test",
        choices=cfg.ALLOWED_CATALOGS,
    )
    parser.add_argument(
        "--flux",
        "-f",
        type=str,
        default=cfg.FLUX_CHOICES[cfg.FALSE_INDEX],
        help="weight the sources with x-ray flux, instead of using the redshift. Possible only with Turin catalog.",
        choices=cfg.FLUX_CHOICES,
    )
    args = parser.parse_args()
    reco = args.reco
    catalog = args.catalog
    flux = args.flux
    if flux == cfg.FLUX_CHOICES[cfg.TRUE_INDEX]:
        flux = True
    elif flux == cfg.FLUX_CHOICES[cfg.FALSE_INDEX]:
        flux = False
    if flux and catalog == cfg.ALLOWED_CATALOGS[cfg.MILLIQUAS_INDEX]:
        raise ValueError(
            f"Possible to use the x-ray fluxes as weighting only with the Turin catalog."
        )
    
    data_path, data_results_path = define_paths()

    filename = None
    url = None
    if catalog == cfg.ALLOWED_CATALOGS[cfg.TURIN_INDEX]:
        filename = cfg.TURIN_FILENAME
        url = cfg.TURIN_URL
    elif catalog == cfg.ALLOWED_CATALOGS[cfg.MILLIQUAS_INDEX]:
        filename = cfg.MILLIQUAS_FILENAME
        url = cfg.MILLIQUAS_URL

    print(f"Checking if '{filename}' is in '{data_path}'...")

    if os.path.isfile(data_path / filename):
        print(f"'{filename}' in '{data_path}', no need to download")

    else:
        print(f"{filename} not found, download {catalog} catalog...")

        r = requests.get(url, allow_redirects=True)
        if catalog == cfg.ALLOWED_CATALOGS[cfg.TURIN_INDEX]:
            catalog_path = data_path / cfg.TURIN_FILENAME
            open(catalog_path, "wb").write(r.content)
        elif catalog == cfg.ALLOWED_CATALOGS[cfg.MILLIQUAS_INDEX]:
            zip_path = data_path / cfg.MILLIQUAS_ZIP
            open(zip_path, "wb").write(r.content)

            print(f"Unzipping {catalog} catalog...")

            with zipfile.ZipFile(zip_path, "r") as zip_ref:
                zip_ref.extractall(data_path)

            print(f"Removing {cfg.MILLIQUAS_ZIP}...")

            os.remove(zip_path)

    names_path = None
    ras_catalog = np.array([])
    decs_catalog = np.array([])
    redshifts_catalog = np.array([])
    names_catalog = np.array([])
    if catalog == cfg.ALLOWED_CATALOGS[cfg.TURIN_INDEX]:
        names_filename = cfg.TURIN_NAMES_FILENAME
        names_path = data_path / names_filename

        print(f"Checking if '{names_filename}' is in '{data_path}'...")

        if os.path.isfile(names_path):
            print(f"'{names_filename}' in '{data_path}', no need to download")
        else:
            url_names = cfg.TURIN_NAMES_URL
            print(f"{names_filename} not found, download from {url_names}...")
            r_names = requests.get(url_names, allow_redirects=True)
            open(names_path, "wb").write(r_names.content)

        if flux:
            xray_filename = cfg.BASS_XRAY_FILENAME
            xray_path = data_path / xray_filename

            print(f"Checking if '{xray_filename}' is in '{data_path}'...")

            if os.path.isfile(xray_path):
                print(f"'{xray_filename}' in '{data_path}', no need to download")
            else:
                url_xray = cfg.BASS_XRAY_URL
                print(f"{xray_filename} not found, download from {url_xray}...")
                r_xray = requests.get(url_xray, allow_redirects=True)
                open(xray_path, "wb").write(r_xray.content)

    print(f"Loading the {catalog} catalog...")

    if catalog == cfg.ALLOWED_CATALOGS[cfg.TURIN_INDEX]:
        dataframe = pd.read_fwf(
            data_path / cfg.TURIN_FILENAME,
            header=cfg.TURIN_HEADER_PRESENT,
            names=cfg.TURIN_HEADER,
        )
        dataframe_names = pd.read_fwf(
            names_path,
            colspecs=cfg.TURIN_NAMES_COLSPECS,
            header=cfg.TURIN_HEADER_PRESENT,
            names=cfg.TURIN_NAMES_HEADER,
        )
        if flux:
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
            xray_catalog = np.array([])
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
            dataframe[cfg.TURIN_RAh], dataframe[cfg.TURIN_RAm], dataframe[cfg.TURIN_RAs]
        ).to_numpy()
        decs_catalog = dms_to_deg(
            dataframe[cfg.TURIN_DEd], dataframe[cfg.TURIN_DEm], dataframe[cfg.TURIN_DEs]
        ).to_numpy()
        redshifts_catalog = dataframe[cfg.TURIN_z].to_numpy()
        names_catalog = dataframe_names[cfg.TURIN_NAME_SOURCE].to_numpy()
        names_catalog[pd.isna(names_catalog)] = dataframe_names[cfg.TURIN_SWIFT][
            pd.isna(names_catalog)
        ]
        names_catalog[pd.isna(names_catalog)] = dataframe_names[cfg.TURIN_WISE][
            pd.isna(names_catalog)
        ]
    elif catalog == cfg.ALLOWED_CATALOGS[cfg.MILLIQUAS_INDEX]:
        dataframe = pd.read_fwf(
            data_path / cfg.MILLIQUAS_FILENAME,
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

    print("Defining the test statistic...")

    effective_area = np.genfromtxt(data_path / cfg.EFFECTIVE_AREA_FILENAME)
    energy_bins = effective_area[:, cfg.EFFECTIVE_AREA_ENERGY_BINS_INDEX]
    effective_area_array = np.array(
        [
            effective_area[:, cfg.EFFECTIVE_AREA_30_90_DEG_INDEX],
            effective_area[:, cfg.EFFECTIVE_AREA_0_30_DEG_INDEX],
            effective_area[:, cfg.EFFECTIVE_AREA_MIN5_0_DEG_INDEX],
            effective_area[:, cfg.EFFECTIVE_AREA_MIN30_MIN5_DEG_INDEX],
            effective_area[:, cfg.EFFECTIVE_AREA_MIN90_MIN30_DEG_INDEX],
        ]
    )

    def energy_factor(bin_index):
        """
        Estimate energy factor for expected number of detected neutrinos
        """
        e_max = energy_bins[bin_index + 1]
        e_min = energy_bins[bin_index]
        factor = (e_max - e_min) / (e_max * e_min)
        return factor

    def area_energy_factor_calculator(a_index):
        """
        Estimate area and energy factor for expected number of detected neutrinos
        """
        eff_area = effective_area_array[a_index]
        factor = 0
        for k in range(len(energy_bins) - 1):
            element = eff_area[k] * energy_factor(k)
            factor += element
        return factor  # units: m^2 GeV^-1

    area_energy_factors = np.array([])
    for i in range(len(effective_area_array)):
        area_energy_factors = np.append(
            area_energy_factors, area_energy_factor_calculator(i)
        )
    hubble_in_s = cfg.HUBBLE_CONSTANT / cfg.MPC_TO_KM
    days = cfg.MJD_04102023 - cfg.MJD_GOLDBRONZE_START
    seconds = (days / cfg.DAYS_IN_YEAR) * cfg.SECONDS_IN_YEAR

    def expected_nu_from_source(z, dec):
        """
        Given the redshift and the declination of a source, determines the total
        number of expected neutrinos from the source
        """
        area_energy_factor = None
        if 90 >= dec > 30:
            area_energy_factor = area_energy_factors[
                cfg.EFFECTIVE_AREA_30_90_DEG_INDEX - 1
            ]
        elif dec <= 30 and dec > 0:
            area_energy_factor = area_energy_factors[
                cfg.EFFECTIVE_AREA_0_30_DEG_INDEX - 1
            ]
        elif dec <= 0 and dec > -5:
            area_energy_factor = area_energy_factors[
                cfg.EFFECTIVE_AREA_MIN5_0_DEG_INDEX - 1
            ]
        elif dec <= -5 and dec > -30:
            area_energy_factor = area_energy_factors[
                cfg.EFFECTIVE_AREA_MIN30_MIN5_DEG_INDEX - 1
            ]
        elif dec <= -30 and dec >= -90:
            area_energy_factor = area_energy_factors[
                cfg.EFFECTIVE_AREA_MIN90_MIN30_DEG_INDEX - 1
            ]
        constant = (
            (hubble_in_s ** 2)
            * seconds
            / (4 * np.pi * (z ** 2) * (cfg.SPEED_OF_LIGHT ** 2))
        )  # m^-2 * s
        expected_nu = constant * cfg.FLUX_NU * (cfg.E0 ** 2) * area_energy_factor
        return expected_nu

    def flux_contribute(z, dec):
        """
        Given the redshift and the declination of a source, determines the contribution
        to the test statistic related to the neutrino flux of the source
        """
        mu = expected_nu_from_source(z, dec)
        contribute = (
            np.log(0.5) + 2 * np.log(mu) - mu
        )  # Here we assume the limit of low fluxes as valid
        return contribute

    def select_effective_area(dec, energy):
        if 90 >= dec > 30:
            effa = effective_area_array[cfg.EFFECTIVE_AREA_30_90_DEG_INDEX - 1]
        elif dec <= 30 and dec > 0:
            effa = effective_area_array[cfg.EFFECTIVE_AREA_0_30_DEG_INDEX - 1]
        elif dec <= 0 and dec > -5:
            effa = effective_area_array[cfg.EFFECTIVE_AREA_MIN5_0_DEG_INDEX - 1]
        elif dec <= -5 and dec > -30:
            effa = effective_area_array[cfg.EFFECTIVE_AREA_MIN30_MIN5_DEG_INDEX - 1]
        elif dec <= -30 and dec >= -90:
            effa = effective_area_array[cfg.EFFECTIVE_AREA_MIN90_MIN30_DEG_INDEX - 1]
        for index in range(len(energy_bins)):
            next_ebin = energy_bins[index + 1]
            if next_ebin >= energy:
                break
        return effa[index]

    if flux:

        def expected_nu_from_source_xray(xray, dec):
            """
            Given the xray flux and the declination of a source, determines the total
            number of expected neutrinos from the source
            """
            area_energy_factor = None
            if 90 >= dec > 30:
                area_energy_factor = area_energy_factors[
                    cfg.EFFECTIVE_AREA_30_90_DEG_INDEX - 1
                ]
            elif dec <= 30 and dec > 0:
                area_energy_factor = area_energy_factors[
                    cfg.EFFECTIVE_AREA_0_30_DEG_INDEX - 1
                ]
            elif dec <= 0 and dec > -5:
                area_energy_factor = area_energy_factors[
                    cfg.EFFECTIVE_AREA_MIN5_0_DEG_INDEX - 1
                ]
            elif dec <= -5 and dec > -30:
                area_energy_factor = area_energy_factors[
                    cfg.EFFECTIVE_AREA_MIN30_MIN5_DEG_INDEX - 1
                ]
            elif dec <= -30 and dec >= -90:
                area_energy_factor = area_energy_factors[
                    cfg.EFFECTIVE_AREA_MIN90_MIN30_DEG_INDEX - 1
                ]
            expected_nu = (
                cfg.CONSTANT_XRAY
                * xray
                * cfg.ERG_TO_GEV
                * (cfg.E0 ** 2)
                * area_energy_factor
            )
            return expected_nu

        def flux_contribute(xray, dec):
            """
            Given the xray flux and the declination of a source, determines the contribution
            to the test statistic related to the neutrino flux of the source
            """
            mu = expected_nu_from_source_xray(xray, dec)
            contribute = (
                np.log(0.5) + 2 * np.log(mu) - mu
            )  # Here we assume the limit of low fluxes as valid
            return contribute

    def unc_contribute(sigma1, sigma2):
        """
        Contribute to the test statistic related only to the uncertainties of the alerts
        """
        return -2 * np.log(sigma1 * sigma2)

    def angular_dist_score(az_true, zen_true, az_pred, zen_pred):
        """
        calculate the MAE of the angular distance between two directions.
        The two vectors are first converted to cartesian unit vectors,
        and then their scalar product is computed, which is equal to
        the cosine of the angle between the two vectors. The inverse
        cosine (arccos) thereof is then the angle between the two input vectors

        Parameters:
        -----------

        az_true : float (or array thereof)
            true azimuth value(s) in radian
        zen_true : float (or array thereof)
            true zenith value(s) in radian
        az_pred : float (or array thereof)
            predicted azimuth value(s) in radian
        zen_pred : float (or array thereof)
            predicted zenith value(s) in radian

        Returns:
        --------

        dist : float
            mean over the angular distance(s) in radian
        """
        if not (
            np.all(np.isfinite(az_true))
            and np.all(np.isfinite(zen_true))
            and np.all(np.isfinite(az_pred))
            and np.all(np.isfinite(zen_pred))
        ):
            raise ValueError("All arguments must be finite")
        # pre-compute all sine and cosine values
        sa1 = np.sin(az_true)
        ca1 = np.cos(az_true)
        sz1 = np.sin(zen_true)
        cz1 = np.cos(zen_true)
        sa2 = np.sin(az_pred)
        ca2 = np.cos(az_pred)
        sz2 = np.sin(zen_pred)
        cz2 = np.cos(zen_pred)
        # scalar product of the two cartesian vectors (x = sz*ca, y = sz*sa, z = cz)
        scalar_prod = sz1 * sz2 * (ca1 * ca2 + sa1 * sa2) + (cz1 * cz2)
        # scalar product of two unit vectors is always between -1 and 1, this is against nummerical instability
        # that might otherwise occure from the finite precision of the sine and cosine functions
        scalar_prod = np.clip(scalar_prod, -1, 1)
        # convert back to an angle (in radian)
        return np.average(np.abs(np.arccos(scalar_prod)))

    def dir_contribute(ra1, dec1, ra2, dec2, raS, decS, sigma1, sigma2):  # in radiants
        """
        Contribute to the test statistic related to
        how near are the alerts to the source. All
        directions must be given in radiants
        """
        phi1 = angular_dist_score(ra1, dec1 + np.pi / 2.0, raS, decS + np.pi / 2.0)
        phi2 = angular_dist_score(ra2, dec2 + np.pi / 2.0, raS, decS + np.pi / 2.0)
        cont = -0.5 * ((phi1 / sigma1) ** 2 + (phi2 / sigma2) ** 2)
        return cont

    def noise_contribute(dec1, dec2, energy1, energy2):
        """
        Contribute to the test statistic related to
        the null hypothesis
        """
        effa1 = select_effective_area(dec1, energy1)
        effa2 = select_effective_area(dec2, energy2)
        return -np.log(np.cos(dec1) * np.cos(dec2) * effa1 * effa2)

    def test_statistic(
        ra1, dec1, sigma1, energy1, ra2, dec2, sigma2, energy2, raS, decS, z_or_xray
    ):
        """
        Test statistic related to two alerts and a source
        """
        c1 = flux_contribute(z_or_xray, decS)
        c2 = unc_contribute(sigma1, sigma2)
        c3 = dir_contribute(ra1, dec1, ra2, dec2, raS, decS, sigma1, sigma2)
        c4 = noise_contribute(dec1, dec2, energy1, energy2)
        ts = c1 + c2 + c3 + c4
        return ts

    print(f"Retrieving the alerts reconstructed with {reco}...")

    RAs = np.array([])
    DECs = np.array([])
    sigmas = np.array([])
    NAMEs = np.array([])
    ENERGIES = np.array([])
    if reco == cfg.ALLOWED_RECONSTRUCTIONS[cfg.SPLINEMPE_INDEX]:
        rev0 = False
        rev1 = False
        has_rev1 = False
        is_data = False
        rev1_names = np.array([])
        splinempe_f = open(data_path / cfg.SPLINEMPE_FILENAME)
        notice_line_index = 0
        for index, line in enumerate(splinempe_f):
            if line == cfg.SPLINEMPE_GCN_START and index > cfg.SPLINEMPE_INDEX_START:
                notice_line_index = 0
                if index < 100:
                    is_data = True
            if line == cfg.SPLINEMPE_COMMENT_START:
                is_data = False
            if is_data:
                if notice_line_index == 2:
                    rev_number = line.split(">")[1].split("<")[0]
                    if rev_number == "1" or rev_number == "2":
                        rev0 = False
                        rev1 = True
                    if rev_number == "0":
                        rev0 = True
                        rev1 = False
                notice_line_index += 1
                if notice_line_index == 4:
                    date = line.split(">")[1].split("<")[0]
                    year = date.split("/")[0]
                    month = date.split("/")[1]
                    day = date.split("/")[2]
                    name = f"IC{year}{month}{day}A"
                    if rev1:
                        if len(rev1_names) > 0:
                            if name in rev1_names:
                                if name == f"IC{year}{month}{day}B":
                                    continue
                                rev1_names[
                                    np.where(rev1_names == name)
                                ] = f"IC{year}{month}{day}B"
                        rev1_names = np.append(rev1_names, name)
                    elif rev0:
                        if len(NAMEs) > 0:
                            if name in NAMEs:
                                NAMEs[
                                    np.where(NAMEs == name)
                                ] = f"IC{year}{month}{day}B"
                        if (
                            name in rev1_names or name in cfg.SPLINEMPE_EXCEPTIONS
                        ) and not name in cfg.SPLINEMPE_BACKGROUND:
                            NAMEs = np.append(NAMEs, name)
                            has_rev1 = True
                        else:
                            has_rev1 = False
                if rev0 and is_data and has_rev1:
                    if notice_line_index == 7:
                        ra = float(line.split(">")[1].split("<")[0])
                        RAs = np.append(RAs, ra)
                    if notice_line_index == 8:
                        de = float(line.split(">")[1].split("<")[0])
                        DECs = np.append(DECs, de)
                    if notice_line_index == 10:
                        err_50 = float(line.split(">")[1].split("<")[0]) / 60
                        sigma = np.deg2rad(err_50) / cfg.RATIO_50_TO_SIGMA
                        sigmas = np.append(sigmas, sigma)
                    if notice_line_index == 11:
                        energy = (
                            float(line.split(">")[1].split("<")[0]) * cfg.TEV_TO_GEV
                        )
                        ENERGIES = np.append(ENERGIES, energy)
    elif reco == cfg.ALLOWED_RECONSTRUCTIONS[cfg.MILLIPEDE_INDEX]:
        alerts_df = pd.read_csv(data_path / cfg.MILLIPEDE_FILENAME)
        RAs = alerts_df[cfg.MILLIPEDE_RA].to_numpy()
        DECs = alerts_df[cfg.MILLIPEDE_DEC].to_numpy()
        RAs_ERR_PLUS = alerts_df[cfg.MILLIPEDE_RA_PLUS].to_numpy()
        DECs_ERR_PLUS = alerts_df[cfg.MILLIPEDE_DEC_PLUS].to_numpy()
        RAs_ERR_MINUS = alerts_df[cfg.MILLIPEDE_RA_MINUS].to_numpy()
        DECs_ERR_MINUS = alerts_df[cfg.MILLIPEDE_DEC_MINUS].to_numpy()
        NAMEs = alerts_df[cfg.MILLIPEDE_IC_NAME].to_numpy()
        ENERGIES = alerts_df[cfg.MILLIPEDE_ENERGY].to_numpy()

        def millipede_area(index):
            """
            Given the index of the alert, returns the angular area
            """
            phi1_deg = RAs[index] - RAs_ERR_MINUS[index]
            phi2_deg = RAs[index] + RAs_ERR_PLUS[index]
            delta1_deg = DECs[index] - DECs_ERR_MINUS[index]
            delta2_deg = DECs[index] + DECs_ERR_PLUS[index]
            phi1_rad = np.deg2rad(phi1_deg)
            phi2_rad = np.deg2rad(phi2_deg)
            delta1_rad = np.deg2rad(delta1_deg)
            delta2_rad = np.deg2rad(delta2_deg)
            A = (phi2_rad - phi1_rad) * (np.sin(delta2_rad) - np.sin(delta1_rad))
            return A

        sigmas = np.array([])
        for i in range(len(alerts_df)):
            area = millipede_area(i)
            sigma = np.sqrt(area / np.pi) / cfg.RATIO_90_TO_SIGMA
            sigmas = np.append(sigmas, sigma)

    print("Estimating background...")

    total_scramblings = 0
    ang_dist_fast_selection = 0
    search_radius = 0
    if reco == cfg.ALLOWED_RECONSTRUCTIONS[cfg.SPLINEMPE_INDEX]:
        ang_dist_fast_selection = cfg.SPLINEMPE_ANG_DIST_FAST_SELECTION
        search_radius = cfg.SPLINEMPE_SEARCH_RADIUS
        if catalog == cfg.ALLOWED_CATALOGS[cfg.TURIN_INDEX]:
            total_scramblings = cfg.TOTAL_SCRAMBLINGS_SPLINEMPE_TURIN
        elif catalog == cfg.ALLOWED_CATALOGS[cfg.MILLIQUAS_INDEX]:
            total_scramblings = cfg.TOTAL_SCRAMBLINGS_SPLINEMPE_MILLIQUAS
    elif reco == cfg.ALLOWED_RECONSTRUCTIONS[cfg.MILLIPEDE_INDEX]:
        ang_dist_fast_selection = cfg.MILLIPEDE_ANG_DIST_FAST_SELECTION
        search_radius = cfg.MILLIPEDE_SEARCH_RADIUS
        if catalog == cfg.ALLOWED_CATALOGS[cfg.TURIN_INDEX]:
            total_scramblings = cfg.TOTAL_SCRAMBLINGS_MILLIPEDE_TURIN
            if flux:
                total_scramblings = cfg.TOTAL_SCRAMBLINGS_MILLIPEDE_TURIN_XRAY
        elif catalog == cfg.ALLOWED_CATALOGS[cfg.MILLIQUAS_INDEX]:
            total_scramblings = cfg.TOTAL_SCRAMBLINGS_MILLIPEDE_MILLIQUAS
    test_statistic_per_scramble = np.array([])
    names_alerts_per_scramble = np.array([])
    names_source_per_scramble = np.array([])
    t0 = time.time()
    for scrambling_number in range(total_scramblings):
        if (scrambling_number + 1) % 100 == 0:
            sys.stdout.write(
                "\r"
                + f"Scrumble nr {scrambling_number + 1:6}"
                + " of "
                + str(total_scramblings)
                + f". Taken {round(time.time() - t0, 1):6}"
                + f" seconds so far. Still {round((total_scramblings - (scrambling_number + 1)) * (time.time() - t0) / (scrambling_number + 1), 1):6}"
                + " seconds remaining."
            )
        rng = np.random.default_rng(seed=scrambling_number)
        random_ras = rng.uniform(0.0, cfg.ROUND_ANGLE, size=len(RAs))
        test_statistic_per_doublet = np.array([])
        names_alerts_per_doublet = np.array([])
        names_source_per_doublet = np.array([])
        for first_alert_index in range(len(RAs)):
            for second_alert_index in range(first_alert_index, len(RAs)):
                if first_alert_index == second_alert_index:
                    continue
                ra1 = random_ras[first_alert_index]
                ra2 = random_ras[second_alert_index]
                dec1 = DECs[first_alert_index]
                dec2 = DECs[second_alert_index]
                higher_ra = max(ra1, ra2)
                smaller_ra = min(ra1, ra2)
                # First fast selection
                ra_distance = min(higher_ra - smaller_ra, smaller_ra - higher_ra + 360)
                dec_distance = np.abs(dec1 - dec2)
                if (
                    ra_distance > ang_dist_fast_selection
                    or dec_distance > ang_dist_fast_selection
                ):
                    continue
                ra1_rad = np.deg2rad(ra1)
                ra2_rad = np.deg2rad(ra2)
                dec1_rad = np.deg2rad(dec1)
                dec2_rad = np.deg2rad(dec2)
                sigma1 = sigmas[first_alert_index]
                sigma2 = sigmas[second_alert_index]
                energy1 = ENERGIES[first_alert_index]
                energy2 = ENERGIES[second_alert_index]
                # Consider the sources nearest to the alert with best angular resolution
                if sigma1 <= sigma2:
                    search_index = first_alert_index
                else:
                    search_index = second_alert_index
                mask_ra = np.logical_and(
                    ras_catalog < random_ras[search_index] + search_radius,
                    ras_catalog > random_ras[search_index] - search_radius,
                )
                mask_dec = np.logical_and(
                    decs_catalog < DECs[search_index] + search_radius,
                    decs_catalog > DECs[search_index] - search_radius,
                )
                mask_sources = np.logical_and(mask_ra, mask_dec)
                test_statistic_per_source = np.array([])
                RAs_sources_nearby = ras_catalog[mask_sources]
                DECs_sources_nearby = decs_catalog[mask_sources]
                if flux:
                    xray_sources_nearby = xray_catalog[mask_sources]
                else:
                    redshifts_sources_nearby = redshifts_catalog[mask_sources]
                names_sources_nearby = names_catalog[mask_sources]
                for source_index in range(len(names_sources_nearby)):
                    ra_rad_source = np.deg2rad(RAs_sources_nearby[source_index])
                    de_rad_source = np.deg2rad(DECs_sources_nearby[source_index])
                    if flux:
                        redshift_or_flux = xray_sources_nearby[source_index]
                    else:
                        redshift_or_flux = redshifts_sources_nearby[source_index]
                    test_statistic_source = test_statistic(
                        ra1_rad,
                        dec1_rad,
                        sigma1,
                        energy1,
                        ra2_rad,
                        dec2_rad,
                        sigma2,
                        energy2,
                        ra_rad_source,
                        de_rad_source,
                        redshift_or_flux,
                    )
                    test_statistic_per_source = np.append(
                        test_statistic_per_source, test_statistic_source
                    )
                if len(test_statistic_per_source) == 0:
                    continue
                index_best_ts = np.argmax(test_statistic_per_source)
                test_statistic_doublet = test_statistic_per_source[index_best_ts]
                name_best_source = names_sources_nearby[index_best_ts]
                test_statistic_per_doublet = np.append(
                    test_statistic_per_doublet, test_statistic_doublet
                )
                names_alerts_per_doublet = np.append(
                    names_alerts_per_doublet,
                    f"{NAMEs[first_alert_index]}, {NAMEs[second_alert_index]}",
                )
                names_source_per_doublet = np.append(
                    names_source_per_doublet, name_best_source
                )
        if len(test_statistic_per_doublet) == 0:
            test_statistic_scramble = cfg.TEST_STATISTIC_EMPTY_SCRAMBLE
            names_alerts_scramble = None
            name_source_scramble = None
        else:
            index_best_doublet = np.argmax(test_statistic_per_doublet)
            test_statistic_scramble = test_statistic_per_doublet[index_best_doublet]
            names_alerts_scramble = names_alerts_per_doublet[index_best_doublet]
            name_source_scramble = names_source_per_doublet[index_best_doublet]
        test_statistic_per_scramble = np.append(
            test_statistic_per_scramble, test_statistic_scramble
        )
        names_alerts_per_scramble = np.append(
            names_alerts_per_scramble, names_alerts_scramble
        )
        names_source_per_scramble = np.append(
            names_source_per_scramble, name_source_scramble
        )

    print("Saving the ts distribution under the background hypothesis...")

    test_statistic_filename = f"{cfg.TEST_STATISTIC_FILENAME}_{reco}_{catalog}"
    if flux:
        test_statistic_filename = f"{test_statistic_filename}_xray"
    else:
        test_statistic_filename = f"{test_statistic_filename}_redshift"
    np.save(data_results_path / test_statistic_filename, test_statistic_per_scramble)

    print(f"\nEstimate ts value for {reco} with {catalog}...")

    test_statistic_per_doublet = np.array([])
    names_alerts_per_doublet = np.array([])
    names_source_per_doublet = np.array([])
    for first_alert_index in range(len(RAs)):
        for second_alert_index in range(first_alert_index, len(RAs)):
            if first_alert_index == second_alert_index:
                continue
            ra1 = RAs[first_alert_index]
            ra2 = RAs[second_alert_index]
            dec1 = DECs[first_alert_index]
            dec2 = DECs[second_alert_index]
            higher_ra = max(ra1, ra2)
            smaller_ra = min(ra1, ra2)
            # First fast selection
            ra_distance = min(higher_ra - smaller_ra, smaller_ra - higher_ra + 360)
            dec_distance = np.abs(dec1 - dec2)
            if (
                ra_distance > ang_dist_fast_selection
                or dec_distance > ang_dist_fast_selection
            ):
                continue
            ra1_rad = np.deg2rad(ra1)
            ra2_rad = np.deg2rad(ra2)
            dec1_rad = np.deg2rad(dec1)
            dec2_rad = np.deg2rad(dec2)
            sigma1 = sigmas[first_alert_index]
            sigma2 = sigmas[second_alert_index]
            energy1 = ENERGIES[first_alert_index]
            energy2 = ENERGIES[second_alert_index]
            # Consider the sources nearest to the alert with best angular resolution
            if sigma1 <= sigma2:
                search_index = first_alert_index
            else:
                search_index = second_alert_index
            mask_ra = np.logical_and(
                ras_catalog < RAs[search_index] + search_radius,
                ras_catalog > RAs[search_index] - search_radius,
            )
            mask_dec = np.logical_and(
                decs_catalog < DECs[search_index] + search_radius,
                decs_catalog > DECs[search_index] - search_radius,
            )
            mask_sources = np.logical_and(mask_ra, mask_dec)
            test_statistic_per_source = np.array([])
            RAs_sources_nearby = ras_catalog[mask_sources]
            DECs_sources_nearby = decs_catalog[mask_sources]
            if flux:
                xray_sources_nearby = xray_catalog[mask_sources]
            else:
                redshifts_sources_nearby = redshifts_catalog[mask_sources]
            names_sources_nearby = names_catalog[mask_sources]
            for source_index in range(len(names_sources_nearby)):
                ra_rad_source = np.deg2rad(RAs_sources_nearby[source_index])
                de_rad_source = np.deg2rad(DECs_sources_nearby[source_index])
                if flux:
                    redshift_or_flux = xray_sources_nearby[source_index]
                else:
                    redshift_or_flux = redshifts_sources_nearby[source_index]
                test_statistic_source = test_statistic(
                    ra1_rad,
                    dec1_rad,
                    sigma1,
                    energy1,
                    ra2_rad,
                    dec2_rad,
                    sigma2,
                    energy2,
                    ra_rad_source,
                    de_rad_source,
                    redshift_or_flux,
                )
                test_statistic_per_source = np.append(
                    test_statistic_per_source, test_statistic_source
                )
            if len(test_statistic_per_source) == 0:
                continue
            index_best_ts = np.argmax(test_statistic_per_source)
            test_statistic_doublet = test_statistic_per_source[index_best_ts]
            name_best_source = names_sources_nearby[index_best_ts]
            test_statistic_per_doublet = np.append(
                test_statistic_per_doublet, test_statistic_doublet
            )
            names_alerts_per_doublet = np.append(
                names_alerts_per_doublet,
                f"{NAMEs[first_alert_index]}, {NAMEs[second_alert_index]}",
            )
            names_source_per_doublet = np.append(
                names_source_per_doublet, name_best_source
            )
    best_test_statistic_index = np.argmax(test_statistic_per_doublet)
    best_test_statistic = test_statistic_per_doublet[best_test_statistic_index]
    best_source = names_source_per_doublet[best_test_statistic_index]
    best_alerts = names_alerts_per_doublet[best_test_statistic_index]

    print(
        f"\nTest statistic for {reco} with {catalog}: {best_test_statistic}\nSource: {best_source}. Doublet: {best_alerts}"
    )

    print("Saving test_statistic result...")

    test_statistic_filename_result = f"{test_statistic_filename}_result"
    np.save(data_results_path / test_statistic_filename_result, best_test_statistic)

    p_value = len(
        test_statistic_per_scramble[test_statistic_per_scramble > best_test_statistic]
    ) / len(test_statistic_per_scramble)

    print(
        f"\n\n*********************\n\n"
        f"p-value: {p_value} = {p_value*100}%"
        f"\n\n*********************\n\n"
    )


if __name__ == "__main__":
    main()
