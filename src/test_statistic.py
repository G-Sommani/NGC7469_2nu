from loading_functions import Loader
import config as cfg
import numpy as np
from scipy.stats import poisson, norm, multinomial  # type: ignore
from scipy.integrate import trapezoid  # type: ignore
from typing import Tuple, List, Union
import catalogs


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


class TestStatistic:

    hubble_in_s = cfg.HUBBLE_CONSTANT / cfg.MPC_TO_KM
    days = cfg.MJD_04102023 - cfg.MJD_GOLDBRONZE_START
    seconds = (days / cfg.DAYS_IN_YEAR) * cfg.SECONDS_IN_YEAR

    def __init__(self, flux: str = cfg.FLUX_CHOICES[cfg.REDSHIFT_INDEX]) -> None:

        print("Defining the test statistic...")

        self.flux = flux
        loader = Loader()
        effective_area = loader.load_effective_area()
        self.energy_bins = effective_area[:, cfg.EFFECTIVE_AREA_ENERGY_BINS_INDEX]
        self.effective_area_array = np.array(
            [
                effective_area[:, cfg.EFFECTIVE_AREA_30_90_DEG_INDEX],
                effective_area[:, cfg.EFFECTIVE_AREA_0_30_DEG_INDEX],
                effective_area[:, cfg.EFFECTIVE_AREA_MIN5_0_DEG_INDEX],
                effective_area[:, cfg.EFFECTIVE_AREA_MIN30_MIN5_DEG_INDEX],
                effective_area[:, cfg.EFFECTIVE_AREA_MIN90_MIN30_DEG_INDEX],
            ]
        )
        self.area_energy_factors = np.array([])
        for i in range(len(self.effective_area_array)):
            self.area_energy_factors = np.append(
                self.area_energy_factors, self.area_energy_factor_calculator(i)
            )
        icecat_decs, icecat_ergs = loader.load_icecat()
        self.bkg_probs_dec_erg = self.estimate_probs_bins_dec_energy(
            icecat_decs, icecat_ergs
        )
        self.expected_nus_per_source = np.array([])
        self.weights_per_source = np.array([])
        self.effa_integral = self.calculate_effa_integral()

    def estimate_probs_bins_dec_energy(
        self,
        icecat_decs: np.ndarray,
        icecat_ergs: np.ndarray,
    ) -> np.ndarray:
        """
        Given the icecat catalog, estimate for each bin of the
        effective area the probability for detecting a neutrino
        in that bin
        """
        dec_bins = cfg.EFFECTIVE_AREA_DEC_BINS
        e_bins = self.energy_bins / 1e3
        n_alerts = len(icecat_decs)
        n_alerts_per_bin = list()
        for e_i in range(len(e_bins) - 1):
            n_alerts_e = list()
            e_bin_down = e_bins[e_i]
            e_bin_up = e_bins[e_i + 1]
            e_bin_mask = np.logical_and(
                icecat_ergs > e_bin_down, icecat_ergs < e_bin_up
            )
            for dec_i in range(len(dec_bins) - 1):
                dec_bin_down = dec_bins[dec_i]
                dec_bin_up = dec_bins[dec_i + 1]
                dec_bin_mask = np.logical_and(
                    icecat_decs > dec_bin_down, icecat_decs < dec_bin_up
                )
                e_dec_bin_mask = np.logical_and(e_bin_mask, dec_bin_mask)
                n_alerts_in_bin = np.sum(e_dec_bin_mask)
                n_alerts_e.append(n_alerts_in_bin)
            n_alerts_per_bin.append(n_alerts_e)
        n_alerts_per_bin = np.array(n_alerts_per_bin)  # type: ignore

        self.n_alerts_per_bin = n_alerts_per_bin

        prob_alerts = list()
        for e_i in range(len(e_bins) - 1):
            prob_alerts_e = list()
            e_bin_down = e_bins[e_i]
            e_bin_up = e_bins[e_i + 1]
            # 'Area' in energy space for an energetically uniform
            # flux of neutrinos
            e_area = (e_bin_up - e_bin_down) / (e_bin_up * e_bin_down)
            e_area *= 1e3
            for dec_i in range(len(dec_bins) - 1):
                dec_bin_down = dec_bins[dec_i]
                dec_bin_up = dec_bins[dec_i + 1]
                # Distribute the neutrinos uniformally over the declination inside the bin
                dec_sin_area = np.sin(np.deg2rad(dec_bin_up)) - np.sin(
                    np.deg2rad(dec_bin_down)
                )
                n_alerts_in_bin = n_alerts_per_bin[e_i][dec_i]
                fraction_alerts = n_alerts_in_bin / n_alerts
                prob_alerts_bin = fraction_alerts / (e_area * dec_sin_area)
                prob_alerts_e.append(prob_alerts_bin)
            prob_alerts.append(prob_alerts_e)

        return np.array(prob_alerts)

    def select_bkg_prob_bin(self, dec: float, energy: float) -> float:
        """
        Given a declination and energy, select the bin in declination and
        energy to tell which is the probability of falling in that bin
        """
        dec_deg = np.rad2deg(dec)
        dec_bins = cfg.EFFECTIVE_AREA_DEC_BINS
        for dec_i in range(len(dec_bins) - 1):
            dec_down = dec_bins[dec_i]
            dec_up = dec_bins[dec_i + 1]
            if np.logical_and(dec_down <= dec_deg, dec_up >= dec_deg):
                dec_index = dec_i
                break
        for e_index in range(len(self.energy_bins)):
            next_ebin = self.energy_bins[e_index + 1]
            if next_ebin >= energy:
                break
        bkg_prob_bin = self.bkg_probs_dec_erg[e_index, dec_index]
        return bkg_prob_bin

    def calculate_effa_integral(self) -> float:
        energy_effa_integrals = list()
        for effa in zip(self.effective_area_array):
            energy_effa_integrals.append(trapezoid(effa, x=self.energy_bins))
        theta_mins = [30, 0, -5, -30, -90]
        theta_maxs = [90, 30, 0, -5, -30]
        factors = list()
        for A, thetam, thetaM in zip(energy_effa_integrals, theta_mins, theta_maxs):
            factors.append(
                A * (np.sin(np.deg2rad(thetaM)) - np.sin(np.deg2rad(thetam)))
            )
        return sum(factors)

    def energy_factor(self, bin_index: int) -> float:
        """
        Estimate energy factor for expected number of detected neutrinos
        """
        e_max = self.energy_bins[bin_index + 1]
        e_min = self.energy_bins[bin_index]
        factor = (e_max - e_min) / (e_max * e_min)
        return factor

    def area_energy_factor_calculator(self, a_index: int) -> float:
        """
        Estimate area and energy factor for expected number of detected neutrinos
        """
        eff_area = self.effective_area_array[a_index]
        factor = 0
        for k in range(len(self.energy_bins) - 1):
            element = eff_area[k] * self.energy_factor(k)
            factor += element
        return factor  # units: m^2 GeV^-1

    def get_area_energy_factor(self, dec: float) -> float:
        dec_deg = np.rad2deg(dec)
        if 90 >= dec_deg > 30:
            area_energy_factor = self.area_energy_factors[
                cfg.EFFECTIVE_AREA_30_90_DEG_INDEX - 1
            ]
        elif dec_deg <= 30 and dec_deg > 0:
            area_energy_factor = self.area_energy_factors[
                cfg.EFFECTIVE_AREA_0_30_DEG_INDEX - 1
            ]
        elif dec_deg <= 0 and dec_deg > -5:
            area_energy_factor = self.area_energy_factors[
                cfg.EFFECTIVE_AREA_MIN5_0_DEG_INDEX - 1
            ]
        elif dec_deg <= -5 and dec_deg > -30:
            area_energy_factor = self.area_energy_factors[
                cfg.EFFECTIVE_AREA_MIN30_MIN5_DEG_INDEX - 1
            ]
        elif dec_deg <= -30 and dec_deg >= -90:
            area_energy_factor = self.area_energy_factors[
                cfg.EFFECTIVE_AREA_MIN90_MIN30_DEG_INDEX - 1
            ]

        return area_energy_factor

    def expected_nu_from_source(self, z_or_xray: float, dec: float) -> float:
        """
        Given the redshift (or the xray flux) and the declination of a source, determines the total
        number of expected neutrinos from the source
        """
        area_energy_factor = self.get_area_energy_factor(dec)
        if self.flux == cfg.FLUX_CHOICES[cfg.XRAY_INDEX]:
            expected_nu = (
                cfg.CONSTANT_XRAY
                * z_or_xray
                * cfg.ERG_TO_GEV
                * (cfg.E0 ** 2)
                * area_energy_factor
            )
        elif self.flux == cfg.FLUX_CHOICES[cfg.REDSHIFT_INDEX]:
            constant = (
                (type(self).hubble_in_s ** 2)
                * type(self).seconds
                / (4 * np.pi * (z_or_xray ** 2) * (cfg.SPEED_OF_LIGHT ** 2))
            )  # m^-2 * s
            expected_nu = constant * cfg.FLUX_NU * (cfg.E0 ** 2) * area_energy_factor
        return expected_nu

    def prob_doublet_from_source(self, z_or_xray: float, dec: float) -> float:
        mu = self.expected_nu_from_source(z_or_xray, dec)
        prob = 1 - (1 + mu) * np.exp(-mu)
        return prob

    def prob_dist_from_source(
        self,
        ra1: float,
        de1: float,
        sigma1: float,
        ra2: float,
        de2: float,
        sigma2: float,
        raS: float,
        deS: float,
    ) -> float:
        const = 1 / (4 * (np.pi ** 2) * (sigma1 ** 2) * (sigma2 ** 2))
        exponent = self.dir_contribute(ra1, de1, ra2, de2, raS, deS, sigma1, sigma2)
        return const * np.exp(exponent)

    def prob_doublet_signal(
        self,
        ra1: float,
        de1: float,
        sigma1: float,
        ra2: float,
        de2: float,
        sigma2: float,
        z_or_xray: float,
        raS: float,
        deS: float,
    ) -> float:
        if self.flux == cfg.FLUX_CHOICES[cfg.NOWEIGHT_INDEX]:
            prob_source = 1.0
        else:
            prob_source = self.prob_doublet_from_source(z_or_xray, deS)
        prob_dist = self.prob_dist_from_source(
            ra1, de1, sigma1, ra2, de2, sigma2, raS, deS
        )
        return prob_source * prob_dist

    def prob_neutrino_background(self, de: float, energy: float) -> float:
        """
        Probabilty for a neutrino uncorrelated with any source
        of being detected at exactly a specific location and
        with a specific energy
        """
        bkg_prob_bin = self.select_bkg_prob_bin(de, energy)
        ra_factor = 1 / (2 * np.pi)
        de_factor = np.cos(de)
        energy_factor = energy ** (-2)
        return ra_factor * de_factor * energy_factor * bkg_prob_bin

    def prob_doublet_background(
        self, de1: float, energy1: float, de2: float, energy2: float
    ) -> float:
        """
        Probability for a doublet uncorrelated with any source of being
        detected at exactly their locations and with their specific
        energies
        """
        bkg_prob1 = self.prob_neutrino_background(de1, energy1)
        bkg_prob2 = self.prob_neutrino_background(de2, energy2)
        return bkg_prob1 * bkg_prob2

    """
    Old version of probability for the background

    def prob_doublet_background(
        self, de1: float, energy1: float, de2: float, energy2: float
    ) -> float:
        effa1 = self.select_effective_area(
            de1, energy1
        ) * energy1 / self.effa_integral
        effa2 = self.select_effective_area(
            de2, energy2
        ) * energy2 / self.effa_integral
        return np.cos(de1) * np.cos(de2) * effa1 * effa2 / (4 * np.pi ** 2)
    """

    def set_weights_catalog(self, catalog: catalogs.Catalog) -> None:
        for index in range(len(catalog.names_catalog)):
            if self.flux == cfg.FLUX_CHOICES[cfg.XRAY_INDEX]:
                prob = self.prob_doublet_from_source(
                    catalog.xray_catalog[index], np.deg2rad(catalog.decs_catalog[index])
                )
            elif self.flux == cfg.FLUX_CHOICES[cfg.REDSHIFT_INDEX]:
                prob = self.prob_doublet_from_source(
                    catalog.redshifts_catalog[index],
                    np.deg2rad(catalog.decs_catalog[index]),
                )
            self.weights_per_source = np.append(self.weights_per_source, prob)
        self.weights_per_source = self.weights_per_source / np.sum(
            self.weights_per_source
        )

    def set_expected_nus_catalog(self, catalog: catalogs.Catalog) -> None:
        for index in range(len(catalog.names_catalog)):
            if self.flux == cfg.FLUX_CHOICES[cfg.XRAY_INDEX]:
                n_nu = self.expected_nu_from_source(
                    catalog.xray_catalog[index], np.deg2rad(catalog.decs_catalog[index])
                )
            elif self.flux == cfg.FLUX_CHOICES[cfg.REDSHIFT_INDEX]:
                n_nu = self.expected_nu_from_source(
                    catalog.redshifts_catalog[index],
                    np.deg2rad(catalog.decs_catalog[index]),
                )
            self.expected_nus_per_source = np.append(self.expected_nus_per_source, n_nu)

    def select_sources_randomly(
        self, random_generator: np.random.Generator, size: int
    ) -> np.ndarray:
        source_indexes = list()
        for i in range(size):
            source_index = np.where(
                multinomial.rvs(
                    1, self.weights_per_source, random_state=random_generator, size=1
                )
                == 1
            )[1][0]
            source_indexes.append(source_index)
        source_indexes = np.asarray(source_indexes)  # type: ignore
        return source_indexes  # type: ignore

    def provide_available_nu_indexes(
        self, dec_source: float, ene_nus: np.ndarray
    ) -> np.ndarray:
        """
        Given the declination of the source, provide the indexes
        of the neutrinos which could be injected without falling in
        an empty bin of the background probability
        """
        dist_source = cfg.MAX_DIST_FROM_SOURCE
        dec_source_down = dec_source - dist_source
        dec_source_up = dec_source + dist_source
        dec_bins = cfg.EFFECTIVE_AREA_DEC_BINS
        for dec_i in range(len(dec_bins) - 1):
            dec_down = dec_bins[dec_i]
            dec_up = dec_bins[dec_i + 1]
            if np.logical_and(dec_down <= dec_source_down, dec_up >= dec_source_down):
                dec_index_down = dec_i
            if np.logical_and(dec_down <= dec_source_up, dec_up >= dec_source_up):
                dec_index_up = dec_i
        dec_indexes = np.linspace(
            dec_index_down, dec_index_up, dec_index_up - dec_index_down + 1, dtype=int
        )
        probs_energy_bins = self.bkg_probs_dec_erg[:, dec_indexes]
        not_empty_energy_bins_indexes = list()
        for energy_prob_index, probs_array in enumerate(probs_energy_bins):
            if np.sum(probs_array == 0.0) == 0:
                not_empty_energy_bins_indexes.append(energy_prob_index)
        not_empty_energy_bins_indexes = np.array(  # type: ignore
            not_empty_energy_bins_indexes
        )
        available_nu_indexes = np.array([])
        for energy_index in not_empty_energy_bins_indexes:
            energy_down = self.energy_bins[energy_index]
            energy_up = self.energy_bins[energy_index + 1]
            indexes_nus_within_energies = np.where(
                np.logical_and(energy_down <= ene_nus, energy_up > ene_nus)
            )[0]
            available_nu_indexes = np.append(
                available_nu_indexes, indexes_nus_within_energies
            )
        return available_nu_indexes

    def select_neutrinos_randomly(
        self,
        dec_source: float,
        ene_nus: np.ndarray,
        random_generator: np.random.Generator,
    ) -> Tuple[int, int]:
        """
        Function to select randomly a doublet making sure that the neutrinos
        will avoid empty bins in the declination-energy space of the
        background probability
        """
        available_nu_indexes = self.provide_available_nu_indexes(dec_source, ene_nus)
        len_available_indexes = len(available_nu_indexes)
        first_index = random_generator.integers(
            low=0, high=len_available_indexes, size=1
        )
        second_index = random_generator.integers(
            low=0, high=len_available_indexes, size=1
        )
        while first_index == second_index:
            second_index = random_generator.integers(
                low=0, high=len_available_indexes, size=1
            )
        first_alert_index = int(available_nu_indexes[first_index][0])
        second_alert_index = int(available_nu_indexes[second_index][0])
        return first_alert_index, second_alert_index

    def select_neutrino_randomly(
        self,
        dec_source: float,
        ene_nus: np.ndarray,
        random_generator: np.random.Generator,
    ) -> int:
        """
        Function to select randomly a neutrino making sure that the neutrinos
        will avoid empty bins in the declination-energy space of the
        background probability
        """
        available_nu_indexes = self.provide_available_nu_indexes(dec_source, ene_nus)
        len_available_indexes = len(available_nu_indexes)
        alert_index = random_generator.integers(
            low=0, high=len_available_indexes, size=1
        )
        return int(available_nu_indexes[alert_index][0])

    def gen_nu_for_sources(
        self, random_state: np.random.Generator
    ) -> Union[List[int], int]:
        return poisson.rvs(mu=self.expected_nus_per_source, random_state=random_state)

    def flux_contribute(self, z_or_xray: float, dec: float) -> float:
        if self.flux == cfg.FLUX_CHOICES[cfg.XRAY_INDEX]:
            """
            Given the redshift (or the xray flux) and the declination of a source, determines the contribution
            to the test statistic related to the neutrino flux of the source
            """
            mu = self.expected_nu_from_source(z_or_xray, dec)
            contribute = (
                np.log(0.5) + 2 * np.log(mu) - mu
            )  # Here we assume the limit of low fluxes as valid

        elif self.flux == cfg.FLUX_CHOICES[cfg.REDSHIFT_INDEX]:
            """
            Given the redshift and the declination of a source, determines the contribution
            to the test statistic related to the neutrino flux of the source
            """
            mu = self.expected_nu_from_source(z_or_xray, dec)
            contribute = (
                np.log(0.5) + 2 * np.log(mu) - mu
            )  # Here we assume the limit of low fluxes as valid
        return contribute

    def select_effective_area(self, dec: float, energy: float) -> float:
        dec_deg = np.rad2deg(dec)
        if 90 >= dec_deg > 30:
            effa = self.effective_area_array[cfg.EFFECTIVE_AREA_30_90_DEG_INDEX - 1]
        elif dec_deg <= 30 and dec > 0:
            effa = self.effective_area_array[cfg.EFFECTIVE_AREA_0_30_DEG_INDEX - 1]
        elif dec_deg <= 0 and dec > -5:
            effa = self.effective_area_array[cfg.EFFECTIVE_AREA_MIN5_0_DEG_INDEX - 1]
        elif dec_deg <= -5 and dec > -30:
            effa = self.effective_area_array[
                cfg.EFFECTIVE_AREA_MIN30_MIN5_DEG_INDEX - 1
            ]
        elif dec_deg <= -30 and dec_deg >= -90:
            effa = self.effective_area_array[
                cfg.EFFECTIVE_AREA_MIN90_MIN30_DEG_INDEX - 1
            ]
        for index in range(len(self.energy_bins)):
            next_ebin = self.energy_bins[index + 1]
            if next_ebin >= energy:
                break
        return effa[index]

    @staticmethod
    def unc_contribute(sigma1: float, sigma2: float) -> float:
        """
        Contribute to the test statistic related only to the uncertainties of the alerts
        """
        return -2 * np.log(sigma1 * sigma2)

    @staticmethod
    def dir_contribute(
        ra1: float,
        dec1: float,
        ra2: float,
        dec2: float,
        raS: float,
        decS: float,
        sigma1: float,
        sigma2: float,
    ) -> float:  # in radiants
        """
        Contribute to the test statistic related to
        how near are the alerts to the source. All
        directions must be given in radiants
        """
        phi1 = angular_dist_score(ra1, dec1 + np.pi / 2.0, raS, decS + np.pi / 2.0)
        phi2 = angular_dist_score(ra2, dec2 + np.pi / 2.0, raS, decS + np.pi / 2.0)
        cont = -0.5 * ((phi1 / sigma1) ** 2 + (phi2 / sigma2) ** 2)
        return cont

    @staticmethod
    def gen_rand_dist(
        raS: float | np.ndarray,
        decS: float | np.ndarray,
        sigma: float | np.ndarray,
        random_state: np.random.Generator,
    ) -> Tuple[float | np.ndarray, float | np.ndarray]:
        nu_ra = norm.rvs(loc=raS, scale=sigma, random_state=random_state)
        nu_dec = norm.rvs(loc=decS, scale=sigma, random_state=random_state)
        return nu_ra, nu_dec

    def noise_contribute(
        self, dec1: float, dec2: float, energy1: float, energy2: float
    ) -> float:
        """
        Contribute to the test statistic related to
        the null hypothesis
        """
        ts1 = np.log(self.prob_neutrino_background(dec1, energy1) * (2 * np.pi))
        ts2 = np.log(self.prob_neutrino_background(dec2, energy2) * (2 * np.pi))
        return -(ts1 + ts2)

    """
    Old noise contribution to the test statistic

    def noise_contribute(
        self, dec1: float, dec2: float, energy1: float, energy2: float
    ) -> float:
        """ """
        Contribute to the test statistic related to
        the null hypothesis
        """ """
        effa1 = self.select_effective_area(dec1, energy1)
        effa2 = self.select_effective_area(dec2, energy2)
        return -np.log(np.cos(dec1) * np.cos(dec2) * effa1 * effa2)
    """

    def test_statistic(
        self,
        ra1: float,
        dec1: float,
        sigma1: float,
        energy1: float,
        ra2: float,
        dec2: float,
        sigma2: float,
        energy2: float,
        raS: float,
        decS: float,
        z_or_xray: float,
    ) -> float:
        """
        Test statistic related to two alerts and a source
        """
        if self.flux == cfg.FLUX_CHOICES[cfg.NOWEIGHT_INDEX]:
            c1 = 0.0
        else:
            c1 = self.flux_contribute(z_or_xray, decS)
        c2 = self.unc_contribute(sigma1, sigma2)
        c3 = self.dir_contribute(ra1, dec1, ra2, dec2, raS, decS, sigma1, sigma2)
        c4 = self.noise_contribute(dec1, dec2, energy1, energy2)
        ts = c1 + c2 + c3 + c4
        return ts
