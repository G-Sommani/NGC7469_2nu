import sys
import numpy as np
import time
import argparse
import config as cfg
from loading_functions import Loader
import catalogs
from test_statistic import TestStatistic
import recos
from typing import Tuple


def estimate_background(
    catalog: catalogs.Catalog,
    reco: recos.Reco,
    test_stat: TestStatistic,
) -> np.ndarray:

    ras_catalog = catalog.ras_catalog
    decs_catalog = catalog.decs_catalog
    redshifts_catalog = catalog.redshifts_catalog
    names_catalog = catalog.names_catalog
    xray_catalog = catalog.xray_catalog
    flux = catalog.xray

    RAs = reco.RAs
    DECs = reco.DECs
    sigmas = reco.sigmas
    NAMEs = reco.NAMEs
    ENERGIES = reco.ENERGIES

    print("Estimating background...")

    total_scramblings = catalog.total_scrambling_possibilities[
        reco.total_scramblings_index
    ]
    ang_dist_fast_selection = reco.ang_dist_fast_selection
    search_radius = reco.search_radius

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
            names_alerts_scramble = ""
            name_source_scramble = ""
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

        return test_statistic_per_scramble


def perform_test(reco_name: str, catalog_name: str, flux: bool) -> None:

    loader = Loader()
    data_results_path = loader.data_results_path
    catalog = catalogs.initiate_catalog(catalog_name, xray=flux)
    loader.load_catalog(catalog)

    test_stat = TestStatistic(flux=flux)

    reco = recos.initiate_reco(reco_name)
    loader.load_reco_data(reco)

    ras_catalog = catalog.ras_catalog
    decs_catalog = catalog.decs_catalog
    redshifts_catalog = catalog.redshifts_catalog
    names_catalog = catalog.names_catalog
    xray_catalog = catalog.xray_catalog
    flux = catalog.xray

    RAs = reco.RAs
    DECs = reco.DECs
    sigmas = reco.sigmas
    NAMEs = reco.NAMEs
    ENERGIES = reco.ENERGIES
    ang_dist_fast_selection = reco.ang_dist_fast_selection
    search_radius = reco.search_radius

    test_statistic_per_scramble = estimate_background(
        catalog, reco, test_stat
    )

    print("Saving the ts distribution under the background hypothesis...")

    test_statistic_filename = (
        f"{cfg.TEST_STATISTIC_FILENAME}_{reco.reco_name}_{catalog.catalog_name}"
    )
    if flux:
        test_statistic_filename = f"{test_statistic_filename}_xray"
    else:
        test_statistic_filename = f"{test_statistic_filename}_redshift"
    np.save(
        loader.data_results_path / test_statistic_filename, test_statistic_per_scramble
    )

    print(f"\nEstimate ts value for {reco_name} with {catalog_name}...")

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
        f"\nTest statistic for {reco_name} with {catalog_name}: {best_test_statistic}\nSource: {best_source}. Doublet: {best_alerts}"
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
    reco_name = args.reco
    catalog_name = args.catalog
    flux = args.flux
    if flux == cfg.FLUX_CHOICES[cfg.TRUE_INDEX]:
        flux = True
    elif flux == cfg.FLUX_CHOICES[cfg.FALSE_INDEX]:
        flux = False

    perform_test(reco_name, catalog_name, flux)


if __name__ == "__main__":
    main()
