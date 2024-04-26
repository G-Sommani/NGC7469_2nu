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

    test_stat = TestStatistic(flux=flux)

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
                    test_statistic_source = test_stat.test_statistic(
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
                test_statistic_source = test_stat.test_statistic(
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
