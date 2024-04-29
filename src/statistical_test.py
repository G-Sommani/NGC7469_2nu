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


def select_sources_nearby(
    search_index: int, ras: np.ndarray, reco: recos.Reco, catalog: catalogs.Catalog
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:

    mask_ra = np.logical_and(
        catalog.ras_catalog < ras[search_index] + reco.search_radius,
        catalog.ras_catalog > ras[search_index] - reco.search_radius,
    )
    mask_dec = np.logical_and(
        catalog.decs_catalog < reco.DECs[search_index] + reco.search_radius,
        catalog.decs_catalog > reco.DECs[search_index] - reco.search_radius,
    )
    mask_sources = np.logical_and(mask_ra, mask_dec)
    RAs_sources_nearby = catalog.ras_catalog[mask_sources]
    DECs_sources_nearby = catalog.decs_catalog[mask_sources]
    if catalog.xray:
        xray_or_redshifts_sources_nearby = catalog.xray_catalog[mask_sources]
    else:
        xray_or_redshifts_sources_nearby = catalog.redshifts_catalog[mask_sources]
    names_sources_nearby = catalog.names_catalog[mask_sources]

    return (
        RAs_sources_nearby,
        DECs_sources_nearby,
        xray_or_redshifts_sources_nearby,
        names_sources_nearby,
    )


def evaluate_doublet(
    ra1_rad: float,
    ra2_rad: float,
    dec1_rad: float,
    dec2_rad: float,
    sigma1: float,
    sigma2: float,
    energy1: float,
    energy2: float,
    names_sources_nearby: np.ndarray,
    RAs_sources_nearby: np.ndarray,
    DECs_sources_nearby: np.ndarray,
    xray_or_redshifts_sources_nearby: np.ndarray,
    test_stat: TestStatistic,
) -> Tuple[float | bool, str]:

    test_statistic_per_source = np.array([])
    for source_index in range(len(names_sources_nearby)):
        ra_rad_source = np.deg2rad(RAs_sources_nearby[source_index])
        de_rad_source = np.deg2rad(DECs_sources_nearby[source_index])
        redshift_or_flux = xray_or_redshifts_sources_nearby[source_index]
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
        return False, ""

    index_best_ts = np.argmax(test_statistic_per_source)
    test_statistic_doublet = test_statistic_per_source[index_best_ts]
    name_best_source = names_sources_nearby[index_best_ts]

    return test_statistic_doublet, name_best_source


def evaluate_dataset(
    ras: np.ndarray,
    reco: recos.Reco,
    catalog: catalogs.Catalog,
    test_stat: TestStatistic,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:

    test_statistic_per_doublet = np.array([])
    names_alerts_per_doublet = np.array([])
    names_source_per_doublet = np.array([])

    for first_alert_index in range(len(ras)):
        for second_alert_index in range(first_alert_index, len(ras)):
            if first_alert_index == second_alert_index:
                continue
            ra1 = ras[first_alert_index]
            ra2 = ras[second_alert_index]
            dec1 = reco.DECs[first_alert_index]
            dec2 = reco.DECs[second_alert_index]
            higher_ra = max(ra1, ra2)
            smaller_ra = min(ra1, ra2)
            # First fast selection
            ra_distance = min(higher_ra - smaller_ra, smaller_ra - higher_ra + 360)
            dec_distance = np.abs(dec1 - dec2)
            if (
                ra_distance > reco.ang_dist_fast_selection
                or dec_distance > reco.ang_dist_fast_selection
            ):
                continue
            ra1_rad = np.deg2rad(ra1)
            ra2_rad = np.deg2rad(ra2)
            dec1_rad = np.deg2rad(dec1)
            dec2_rad = np.deg2rad(dec2)
            sigma1 = reco.sigmas[first_alert_index]
            sigma2 = reco.sigmas[second_alert_index]
            energy1 = reco.ENERGIES[first_alert_index]
            energy2 = reco.ENERGIES[second_alert_index]

            # Consider the sources nearest to the alert with best angular resolution
            if sigma1 <= sigma2:
                search_index = first_alert_index
            else:
                search_index = second_alert_index

            (
                RAs_sources_nearby,
                DECs_sources_nearby,
                xray_or_redshifts_sources_nearby,
                names_sources_nearby,
            ) = select_sources_nearby(search_index, ras, reco, catalog)

            test_statistic_doublet, name_best_source = evaluate_doublet(
                ra1_rad,
                ra2_rad,
                dec1_rad,
                dec2_rad,
                sigma1,
                sigma2,
                energy1,
                energy2,
                names_sources_nearby,
                RAs_sources_nearby,
                DECs_sources_nearby,
                xray_or_redshifts_sources_nearby,
                test_stat,
            )

            if not test_statistic_doublet:
                continue

            test_statistic_per_doublet = np.append(
                test_statistic_per_doublet, test_statistic_doublet
            )
            names_alerts_per_doublet = np.append(
                names_alerts_per_doublet,
                f"{reco.NAMEs[first_alert_index]}, {reco.NAMEs[second_alert_index]}",
            )
            names_source_per_doublet = np.append(
                names_source_per_doublet, name_best_source
            )

    return (
        test_statistic_per_doublet,
        names_alerts_per_doublet,
        names_source_per_doublet,
    )


def estimate_background(
    catalog: catalogs.Catalog,
    reco: recos.Reco,
    test_stat: TestStatistic,
    hypo: bool,
) -> np.ndarray:

    RAs = reco.RAs

    print("Estimating background...")

    total_scramblings = catalog.total_scrambling_possibilities[
        reco.total_scramblings_index
    ]
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
        if hypo:
            N_nus_sources = np.array([])
            indexes_emitting_sources = np.array([])
            for s_i in range(len(catalog.redshifts_catalog)):
                N_nu = test_stat.gen_nu_from_source(
                    catalog.redshifts_catalog[s_i], catalog.decs_catalog[s_i], rng
                )
                N_nus_sources = np.append(N_nus_sources, N_nu)
                if N_nu != 0:
                    for k in range(N_nu):
                        indexes_emitting_sources = np.append(indexes_emitting_sources, int(s_i))
            number_source_nus = int(sum(N_nus_sources))
            if number_source_nus != 0:
                remove_neutrinos_indexes = rng.integers(low=0, high=len(reco.NAMEs), size = number_source_nus)
                for i, index in enumerate(remove_neutrinos_indexes):
                    print(i, index, indexes_emitting_sources, indexes_emitting_sources[i])
                    nu_ra_rad, nu_dec_rad = test_stat.gen_rand_dist(
                        np.deg2rad(catalog.ras_catalog[int(indexes_emitting_sources[i])]),
                        np.deg2rad(catalog.decs_catalog[int(indexes_emitting_sources[i])]),
                        reco.sigmas[index],
                        random_state=rng
                    )
                    print(
                        catalog.ras_catalog[int(indexes_emitting_sources[i])], 
                        catalog.decs_catalog[int(indexes_emitting_sources[i])], 
                        np.rad2deg(nu_ra_rad), np.rad2deg(nu_dec_rad)
                    )

        (
            test_statistic_per_doublet,
            names_alerts_per_doublet,
            names_source_per_doublet,
        ) = evaluate_dataset(random_ras, reco, catalog, test_stat)

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


def perform_test(reco_name: str, catalog_name: str, flux: bool, hypo: bool) -> None:

    loader = Loader()
    data_results_path = loader.data_results_path
    catalog = catalogs.initiate_catalog(catalog_name, xray=flux)
    loader.load_catalog(catalog)

    test_stat = TestStatistic(flux=flux)

    reco = recos.initiate_reco(reco_name)
    loader.load_reco_data(reco)

    flux = catalog.xray

    RAs = reco.RAs

    test_statistic_per_scramble = estimate_background(catalog, reco, test_stat, hypo)

    if hypo:
        print("Saving the ts distribution under the alternative hypothesis...")
    else:
        print("Saving the ts distribution under the background hypothesis...")

    test_statistic_filename = (
        f"{cfg.TEST_STATISTIC_FILENAME}_{reco.reco_name}_{catalog.catalog_name}"
    )
    if flux:
        test_statistic_filename = f"{test_statistic_filename}_xray"
    else:
        test_statistic_filename = f"{test_statistic_filename}_redshift"
    if hypo:
        test_statistic_filename = f"{test_statistic_filename}_alternative_hypothesis"
    else:
        test_statistic_filename = f"{test_statistic_filename}_background_hypothesis"
    np.save(
        loader.data_results_path / test_statistic_filename, test_statistic_per_scramble
    )

    print(f"\nEstimate ts value for {reco_name} with {catalog_name}...")

    (
        test_statistic_per_doublet,
        names_alerts_per_doublet,
        names_source_per_doublet,
    ) = evaluate_dataset(RAs, reco, catalog, test_stat)

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
        default=cfg.DEFAULT_RECO,
        help="reconstruction to use for the neutrino events",
        choices=cfg.ALLOWED_RECONSTRUCTIONS,
    )
    parser.add_argument(
        "--catalog",
        "-c",
        type=str,
        default=cfg.DEFAULT_CATALOG,
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
    parser.add_argument(
        "--alternative_hypothesis",
        "-a",
        type=str,
        default=cfg.HYPO_CHOICES[cfg.FALSE_INDEX],
        help="Generate the test statistic distribution using the alternative hypothesis.",
        choices=cfg.HYPO_CHOICES,
    )
    args = parser.parse_args()
    reco_name = args.reco
    catalog_name = args.catalog
    flux = args.flux
    if flux == cfg.FLUX_CHOICES[cfg.TRUE_INDEX]:
        flux = True
    elif flux == cfg.FLUX_CHOICES[cfg.FALSE_INDEX]:
        flux = False
    hypo = args.alternative_hypothesis
    if hypo == cfg.HYPO_CHOICES[cfg.TRUE_INDEX]:
        hypo = True
    elif hypo == cfg.HYPO_CHOICES[cfg.FALSE_INDEX]:
        hypo = False

    perform_test(reco_name, catalog_name, flux, hypo)


if __name__ == "__main__":
    main()
