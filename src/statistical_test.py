import sys
import numpy as np
import time
import argparse
import config as cfg
from loading_functions import Loader
import catalogs
from test_statistic import TestStatistic
import recos
from typing import Tuple, Union
import copy


def select_sources_nearby(
    search_index: int,
    ras: np.ndarray,
    decs: np.ndarray,
    reco: recos.Reco,
    catalog: catalogs.Catalog,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:

    mask_ra = np.logical_and(
        catalog.ras_catalog < ras[search_index] + reco.search_radius,
        catalog.ras_catalog > ras[search_index] - reco.search_radius,
    )
    mask_dec = np.logical_and(
        catalog.decs_catalog < decs[search_index] + reco.search_radius,
        catalog.decs_catalog > decs[search_index] - reco.search_radius,
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

    test_statistic_per_source = list()
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
        test_statistic_per_source.append(test_statistic_source)
    test_statistic_per_source = np.asarray(test_statistic_per_source)  # type: ignore
    if len(test_statistic_per_source) == 0:
        return False, ""

    index_best_ts = np.argmax(test_statistic_per_source)
    test_statistic_doublet = test_statistic_per_source[index_best_ts]
    name_best_source = names_sources_nearby[index_best_ts]

    return test_statistic_doublet, name_best_source


def evaluate_dataset(
    ras: np.ndarray,
    decs: np.ndarray,
    reco: recos.Reco,
    catalog: catalogs.Catalog,
    test_stat: TestStatistic,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:

    test_statistic_per_doublet = list()
    names_alerts_per_doublet = list()
    names_source_per_doublet = list()

    for first_alert_index in range(len(ras)):
        for second_alert_index in range(first_alert_index, len(ras)):
            if first_alert_index == second_alert_index:
                continue
            ra1 = ras[first_alert_index]
            ra2 = ras[second_alert_index]
            dec1 = decs[first_alert_index]
            dec2 = decs[second_alert_index]
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
            ) = select_sources_nearby(search_index, ras, decs, reco, catalog)

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

            test_statistic_per_doublet.append(test_statistic_doublet)
            names_alerts_per_doublet.append(
                f"{reco.NAMEs[first_alert_index]}, {reco.NAMEs[second_alert_index]}",
            )
            names_source_per_doublet.append(name_best_source)
    test_statistic_per_doublet = np.asarray(test_statistic_per_doublet)  # type: ignore
    names_alerts_per_doublet = np.asarray(names_alerts_per_doublet)  # type: ignore
    names_source_per_doublet = np.asarray(names_source_per_doublet)  # type: ignore

    return (  # type: ignore
        test_statistic_per_doublet,
        names_alerts_per_doublet,
        names_source_per_doublet,
    )


def estimate_background(
    catalog: catalogs.Catalog,
    reco: recos.Reco,
    test_stat: TestStatistic,
    hypo: str,
    loader: Loader,
    dec_jitter: float | None,
) -> np.ndarray:

    RAs = reco.RAs

    print("Estimating background...")

    total_scramblings = catalog.total_scrambling_possibilities[
        reco.total_scramblings_index
    ]
    test_statistic_per_scramble = list()
    names_alerts_per_scramble = list()
    names_source_per_scramble = list()
    alternative_probs = list()
    background_probs = list()
    t0 = time.time()
    if hypo == cfg.HYPO_CHOICES[cfg.DOUBLET_INDEX]:
        rng = np.random.default_rng(seed=total_scramblings)
        first_alert_indexes = rng.integers(
            low=0, high=len(reco.NAMEs), size=total_scramblings
        )
        second_alert_indexes = rng.integers(
            low=0, high=len(reco.NAMEs), size=total_scramblings
        )
        inj_source_indexes = test_stat.select_sources_randomly(
            random_generator=rng, size=total_scramblings
        )

        ras_sources = catalog.ras_catalog[inj_source_indexes]
        decs_sources = catalog.decs_catalog[inj_source_indexes]
        if test_stat.flux:
            redshifts_or_xray_sources = catalog.xray_catalog[inj_source_indexes]
        else:
            redshifts_or_xray_sources = catalog.redshifts_catalog[inj_source_indexes]
        sigmas_first = reco.sigmas[first_alert_indexes]
        sigmas_second = reco.sigmas[second_alert_indexes]
        energies_first = reco.ENERGIES[first_alert_indexes]
        energies_second = reco.ENERGIES[second_alert_indexes]
        names_first = reco.NAMEs[first_alert_indexes]
        names_second = reco.NAMEs[second_alert_indexes]
        names_sources = catalog.names_catalog[inj_source_indexes]

        nu1_ras_rad, nu1_decs_rad = test_stat.gen_rand_dist(
            np.deg2rad(ras_sources),
            np.deg2rad(decs_sources),
            sigmas_first,
            random_state=rng,
        )
        nu2_ras_rad, nu2_decs_rad = test_stat.gen_rand_dist(
            np.deg2rad(ras_sources),
            np.deg2rad(decs_sources),
            sigmas_second,
            random_state=rng,
        )
        for i in range(total_scramblings):

            if (i + 1) % 1000 == 0:
                sys.stdout.write(
                    "\r"
                    + f"Scrumble nr {i + 1:6}"
                    + " of "
                    + str(total_scramblings)
                    + f". Taken {round(time.time() - t0, 1):6}"
                    + f" seconds so far. Still {round((total_scramblings - (i + 1)) * (time.time() - t0) / (i + 1), 1):6}"
                    + " seconds remaining."
                )

            ra = ras_sources[i]
            de = decs_sources[i]
            redshift = redshifts_or_xray_sources[i]
            sigma_first = sigmas_first[i]
            sigma_second = sigmas_second[i]
            energy_first = energies_first[i]
            energy_second = energies_second[i]
            nu1_ra_rad = nu1_ras_rad[i]  # type: ignore
            nu1_dec_rad = nu1_decs_rad[i]  # type: ignore
            nu2_ra_rad = nu2_ras_rad[i]  # type: ignore
            nu2_dec_rad = nu2_decs_rad[i]  # type: ignore
            name_first = names_first[i]
            name_second = names_second[i]
            name_source = names_sources[i]

            doublet_ts = test_stat.test_statistic(
                nu1_ra_rad,
                nu1_dec_rad,
                sigma_first,
                energy_first,
                nu2_ra_rad,
                nu2_dec_rad,
                sigma_second,
                energy_second,
                np.deg2rad(ra),
                np.deg2rad(de),
                redshift,
            )
            alternative_prob = test_stat.prob_doublet_signal(
                nu1_ra_rad,
                nu1_dec_rad,
                sigma_first,
                nu2_ra_rad,
                nu2_dec_rad,
                sigma_second,
                redshift,
                np.deg2rad(ra),
                np.deg2rad(de),
            )
            background_prob = test_stat.prob_doublet_background(
                nu1_dec_rad, energy_first, nu2_dec_rad, energy_second
            )
            alternative_probs.append(alternative_prob)
            background_probs.append(background_prob)
            test_statistic_per_scramble.append(doublet_ts)
            names_alerts_per_scramble.append(f"{name_first}, {name_second}")
            names_source_per_scramble.append(name_source)
        alternative_probs = np.asarray(alternative_probs)  # type: ignore
        background_probs = np.asarray(background_probs)  # type: ignore
        test_statistic_per_scramble = np.asarray(test_statistic_per_scramble)  # type: ignore
        if test_stat.flux:
            weight_name = "xray"
        else:
            weight_name = "redshift"
        np.save(
            loader.data_results_path
            / f"{weight_name}_sources_{catalog.catalog_name}_{hypo}_hypothesis",
            redshifts_or_xray_sources,
        )
        np.save(
            loader.data_results_path
            / f"decs_sources_{weight_name}_{catalog.catalog_name}_{hypo}_hypothesis",
            decs_sources,
        )
        np.save(
            loader.data_results_path
            / f"background_probs_{reco.reco_name}_{catalog.catalog_name}_{weight_name}_{hypo}_hypothesis",
            background_probs,
        )
        np.save(
            loader.data_results_path
            / f"alternative_probs_{reco.reco_name}_{catalog.catalog_name}_{weight_name}_{hypo}_hypothesis",
            alternative_probs,
        )
        return test_statistic_per_scramble  # type: ignore

    rng = np.random.default_rng(seed=total_scramblings)
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
        modified_decs = copy.copy(reco.DECs)
        if dec_jitter:
            for index in range(len(modified_decs)):
                dec = reco.DECs[index]
                ene = reco.ENERGIES[index]
                dec_jit = rng.uniform(-dec_jitter, dec_jitter, size=1)[0]
                new_dec = dec + dec_jit
                while test_stat.select_bkg_prob_bin(np.deg2rad(new_dec), ene) == 0.0:
                    dec_jit = rng.uniform(-dec_jitter, dec_jitter, size=len(reco.DECs))[
                        0
                    ]
                    new_dec = dec + dec_jit
                modified_decs[index] = new_dec

        if hypo == cfg.HYPO_CHOICES[cfg.POPULATION_INDEX]:
            N_nus_sources: Union[list[int], int, np.ndarray] = list()
            indexes_emitting_sources: Union[list[int], int, np.ndarray] = list()
            N_nus_sources = test_stat.gen_nu_for_sources(rng)
            N_nus_sources = np.asarray(N_nus_sources)
            indexes_emitting_sources = np.where(N_nus_sources != 0)[0]
            number_source_nus = int(sum(N_nus_sources))
            if number_source_nus != 0:
                remove_neutrinos_indexes = rng.integers(
                    low=0, high=len(reco.NAMEs), size=number_source_nus
                )
                remove_index = 0
                for i, source_index in enumerate(indexes_emitting_sources):
                    source_index = int(source_index)
                    for n_nu in range(N_nus_sources[source_index]):
                        nu_index = remove_neutrinos_indexes[remove_index]
                        nu_ra_rad, nu_dec_rad = test_stat.gen_rand_dist(
                            np.deg2rad(catalog.ras_catalog[source_index]),
                            np.deg2rad(catalog.decs_catalog[source_index]),
                            reco.sigmas[nu_index],
                            random_state=rng,
                        )
                        random_ras[nu_index] = np.rad2deg(nu_ra_rad)
                        modified_decs[nu_index] = np.rad2deg(nu_dec_rad)
                        remove_index += 1

        if hypo == cfg.HYPO_CHOICES[cfg.DOUBLET_INJ_INDEX]:
            inj_source_index = test_stat.select_sources_randomly(
                random_generator=rng, size=1
            )
            ra_source = catalog.ras_catalog[inj_source_index]
            dec_source = catalog.decs_catalog[inj_source_index]
            (
                first_alert_index,
                second_alert_index,
            ) = test_stat.select_neutrinos_randomly(
                dec_source,
                reco.ENERGIES,
                random_generator=rng,
            )
            sigma_first = reco.sigmas[first_alert_index]
            sigma_second = reco.sigmas[second_alert_index]
            nu1_ra_rad, nu1_dec_rad = test_stat.gen_rand_dist(
                np.deg2rad(ra_source),
                np.deg2rad(dec_source),
                sigma_first,
                random_state=rng,
            )
            nu2_ra_rad, nu2_dec_rad = test_stat.gen_rand_dist(
                np.deg2rad(ra_source),
                np.deg2rad(dec_source),
                sigma_second,
                random_state=rng,
            )
            random_ras[first_alert_index] = np.rad2deg(nu1_ra_rad)
            random_ras[second_alert_index] = np.rad2deg(nu2_ra_rad)
            modified_decs[first_alert_index] = np.rad2deg(nu1_dec_rad)
            modified_decs[second_alert_index] = np.rad2deg(nu2_dec_rad)

        if hypo == cfg.HYPO_CHOICES[cfg.SINGLET_INJ_INDEX]:
            inj_source_index = test_stat.select_sources_randomly(
                random_generator=rng, size=1
            )
            ra_source = catalog.ras_catalog[inj_source_index]
            dec_source = catalog.decs_catalog[inj_source_index]
            alert_index = test_stat.select_neutrino_randomly(
                dec_source,
                reco.ENERGIES,
                random_generator=rng,
            )
            sigma = reco.sigmas[alert_index]
            nu_ra_rad, nu_dec_rad = test_stat.gen_rand_dist(
                np.deg2rad(ra_source),
                np.deg2rad(dec_source),
                sigma,
                random_state=rng,
            )
            random_ras[alert_index] = np.rad2deg(nu_ra_rad)
            modified_decs[alert_index] = np.rad2deg(nu_dec_rad)

        (
            test_statistic_per_doublet,
            names_alerts_per_doublet,
            names_source_per_doublet,
        ) = evaluate_dataset(random_ras, modified_decs, reco, catalog, test_stat)

        if len(test_statistic_per_doublet) == 0:
            test_statistic_scramble = cfg.TEST_STATISTIC_EMPTY_SCRAMBLE
            names_alerts_scramble = ""
            name_source_scramble = ""
        else:
            index_best_doublet = np.argmax(test_statistic_per_doublet)
            test_statistic_scramble = test_statistic_per_doublet[index_best_doublet]
            names_alerts_scramble = names_alerts_per_doublet[index_best_doublet]
            name_source_scramble = names_source_per_doublet[index_best_doublet]
            name_first, name_second = names_alerts_scramble.split(", ")
            i_1 = np.where(reco.NAMEs == name_first)[0][0]
            i_2 = np.where(reco.NAMEs == name_second)[0][0]
            i_S = np.where(catalog.names_catalog == name_source_scramble)[0][0]
            p_S = test_stat.prob_doublet_signal(
                np.deg2rad(random_ras[i_1]),
                np.deg2rad(reco.DECs[i_1]),
                reco.sigmas[i_1],
                np.deg2rad(random_ras[i_2]),
                np.deg2rad(reco.DECs[i_2]),
                reco.sigmas[i_2],
                catalog.redshifts_catalog[i_S],
                np.deg2rad(catalog.ras_catalog[i_S]),
                np.deg2rad(catalog.decs_catalog[i_S]),
            )
            p_B = test_stat.prob_doublet_background(
                np.deg2rad(reco.DECs[i_1]),
                reco.ENERGIES[i_1],
                np.deg2rad(reco.DECs[i_2]),
                reco.ENERGIES[i_2],
            )
            background_probs.append(p_B)
            alternative_probs.append(p_S)
        test_statistic_per_scramble.append(test_statistic_scramble)
        names_alerts_per_scramble.append(names_alerts_scramble)
        names_source_per_scramble.append(name_source_scramble)

    test_statistic_per_scramble = np.asarray(test_statistic_per_scramble)  # type: ignore
    names_source_per_scramble = np.asarray(names_source_per_scramble)  # type: ignore
    names_alerts_per_scramble = np.asarray(names_alerts_per_scramble)  # type: ignore
    background_probs = np.asarray(background_probs)  # type: ignore
    alternative_probs = np.asarray(alternative_probs)  # type: ignore
    if test_stat.flux:
        weight_name = "xray"
    else:
        weight_name = "redshift"
    if dec_jitter:
        jitter_name = "_dec_jitter"
    else:
        jitter_name = ""
    np.save(
        loader.data_results_path
        / f"background_probs_{reco.reco_name}_{catalog.catalog_name}_{weight_name}_{jitter_name}_{hypo}_hypothesis",
        background_probs,
    )
    np.save(
        loader.data_results_path
        / f"alternative_probs_{reco.reco_name}_{catalog.catalog_name}_{weight_name}{jitter_name}_{hypo}_hypothesis",
        alternative_probs,
    )

    return test_statistic_per_scramble  # type: ignore


def perform_test(
    reco_name: str, catalog_name: str, flux: str, hypo: str, dec_jitter: float | None
) -> None:

    loader = Loader()
    data_results_path = loader.data_results_path
    if flux == cfg.FLUX_CHOICES[cfg.XRAY_INDEX]:
        xray = True
        noweight = False
    elif flux == cfg.FLUX_CHOICES[cfg.NOWEIGHT_INDEX]:
        xray = False
        noweight = True
    else:
        xray = False
        noweight = False
    print(f"The weight choice is {flux} and noweight is {noweight}")
    catalog = catalogs.initiate_catalog(catalog_name, xray=xray, noweight=noweight)
    loader.load_catalog(catalog)

    test_stat = TestStatistic(flux=flux)
    if hypo is not cfg.HYPO_CHOICES[cfg.BACKGROUND_INDEX]:
        test_stat.set_expected_nus_catalog(catalog)
        test_stat.set_weights_catalog(catalog)

    reco = recos.initiate_reco(reco_name)
    loader.load_reco_data(reco)

    RAs = reco.RAs
    DECs = reco.DECs

    test_statistic_per_scramble = estimate_background(
        catalog, reco, test_stat, hypo, loader, dec_jitter
    )

    if hypo:
        print("Saving the ts distribution under the alternative hypothesis...")
    else:
        print("Saving the ts distribution under the background hypothesis...")

    test_statistic_filename = (
        f"{cfg.TEST_STATISTIC_FILENAME}_{reco.reco_name}_{catalog.catalog_name}"
    )
    test_statistic_filename = f"{test_statistic_filename}_{flux}"
    if dec_jitter:
        test_statistic_filename = f"{test_statistic_filename}_dec_jitter_{dec_jitter}"

    test_statistic_filename = f"{test_statistic_filename}_{hypo}_hypothesis"
    np.save(
        loader.data_results_path / test_statistic_filename, test_statistic_per_scramble
    )

    print(f"\nEstimate ts value for {reco_name} with {catalog_name}...")

    (
        test_statistic_per_doublet,
        names_alerts_per_doublet,
        names_source_per_doublet,
    ) = evaluate_dataset(RAs, DECs, reco, catalog, test_stat)

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
        default=cfg.FLUX_CHOICES[cfg.REDSHIFT_INDEX],
        help="weight the sources with x-ray flux, instead of using the redshift. Possible only with Turin catalog.",
        choices=cfg.FLUX_CHOICES,
    )
    parser.add_argument(
        "--alternative_hypothesis",
        "-a",
        type=str,
        default=cfg.HYPO_CHOICES[cfg.BACKGROUND_INDEX],
        help="Generate the test statistic distribution using the alternative hypothesis.",
        choices=cfg.HYPO_CHOICES,
    )
    parser.add_argument(
        "--declination-jitter",
        "-d",
        type=float,
        default=cfg.DEC_JITTER_DEFAULT,
        help="Apply a declination jitter to the neutrino dataset.",
    )
    args = parser.parse_args()
    reco_name = args.reco
    catalog_name = args.catalog
    flux = args.flux
    dec_jitter = args.declination_jitter
    hypo = args.alternative_hypothesis

    print(
        f"\n\nPerform test with the {catalog_name} catalog,\nthe {reco_name} reconstruction,\nunder the {hypo} hypothesis."
    )
    print(f"Use {flux} as weight\n\n")
    if dec_jitter:
        print(f"Apply a jitter to the declination of {dec_jitter} degrees\n\n")
    else:
        print("No jitter to the declination\n\n")

    perform_test(reco_name, catalog_name, flux, hypo, dec_jitter)


if __name__ == "__main__":
    main()
