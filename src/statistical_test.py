import sys
import pandas as pd  # type: ignore
import numpy as np
import time
import argparse
import config as cfg
from loading_functions import Loader
import catalogs
from test_statistic import TestStatistic


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
    catalog_name = args.catalog
    flux = args.flux
    if flux == cfg.FLUX_CHOICES[cfg.TRUE_INDEX]:
        flux = True
    elif flux == cfg.FLUX_CHOICES[cfg.FALSE_INDEX]:
        flux = False

    loader = Loader()
    data_path = loader.data_path
    data_results_path = loader.data_results_path
    catalog = catalogs.initiate_catalog(catalog_name, xray=flux)
    loader.load_catalog(catalog)

    ras_catalog = catalog.ras_catalog
    decs_catalog = catalog.decs_catalog
    redshifts_catalog = catalog.redshifts_catalog
    names_catalog = catalog.names_catalog
    xray_catalog = catalog.xray_catalog

    test_stat = TestStatistic()
    energy_bins = test_stat.energy_bins
    effective_area_array = test_stat.effective_area_array
    area_energy_factor_calculator = test_stat.area_energy_factor_calculator

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
        if catalog_name == cfg.ALLOWED_CATALOGS[cfg.TURIN_INDEX]:
            total_scramblings = cfg.TOTAL_SCRAMBLINGS_SPLINEMPE_TURIN
        elif catalog_name == cfg.ALLOWED_CATALOGS[cfg.MILLIQUAS_INDEX]:
            total_scramblings = cfg.TOTAL_SCRAMBLINGS_SPLINEMPE_MILLIQUAS
    elif reco == cfg.ALLOWED_RECONSTRUCTIONS[cfg.MILLIPEDE_INDEX]:
        ang_dist_fast_selection = cfg.MILLIPEDE_ANG_DIST_FAST_SELECTION
        search_radius = cfg.MILLIPEDE_SEARCH_RADIUS
        if catalog_name == cfg.ALLOWED_CATALOGS[cfg.TURIN_INDEX]:
            total_scramblings = cfg.TOTAL_SCRAMBLINGS_MILLIPEDE_TURIN
            if flux:
                total_scramblings = cfg.TOTAL_SCRAMBLINGS_MILLIPEDE_TURIN_XRAY
        elif catalog_name == cfg.ALLOWED_CATALOGS[cfg.MILLIQUAS_INDEX]:
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

    test_statistic_filename = f"{cfg.TEST_STATISTIC_FILENAME}_{reco}_{catalog_name}"
    if flux:
        test_statistic_filename = f"{test_statistic_filename}_xray"
    else:
        test_statistic_filename = f"{test_statistic_filename}_redshift"
    np.save(data_results_path / test_statistic_filename, test_statistic_per_scramble)

    print(f"\nEstimate ts value for {reco} with {catalog_name}...")

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
        f"\nTest statistic for {reco} with {catalog_name}: {best_test_statistic}\nSource: {best_source}. Doublet: {best_alerts}"
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
