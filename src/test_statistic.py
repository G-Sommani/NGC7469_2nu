from loading_functions import Loader
import config as cfg
import numpy as np
from scipy.stats import poisson, norm, multinomial  # type: ignore
from scipy.integrate import trapezoid  # type: ignore
from typing import Tuple
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

    def __init__(self, flux: bool = False) -> None:

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
        self.weights_per_source = np.array([])
        self.effa_integral = self.calculate_effa_integral()

    def calculate_effa_integral(self) -> float:
        energy_effa_integrals = list()
        for effa in self.effective_area_array:
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

        if 90 >= dec > 30:
            area_energy_factor = self.area_energy_factors[
                cfg.EFFECTIVE_AREA_30_90_DEG_INDEX - 1
            ]
        elif dec <= 30 and dec > 0:
            area_energy_factor = self.area_energy_factors[
                cfg.EFFECTIVE_AREA_0_30_DEG_INDEX - 1
            ]
        elif dec <= 0 and dec > -5:
            area_energy_factor = self.area_energy_factors[
                cfg.EFFECTIVE_AREA_MIN5_0_DEG_INDEX - 1
            ]
        elif dec <= -5 and dec > -30:
            area_energy_factor = self.area_energy_factors[
                cfg.EFFECTIVE_AREA_MIN30_MIN5_DEG_INDEX - 1
            ]
        elif dec <= -30 and dec >= -90:
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
        if self.flux:
            expected_nu = (
                cfg.CONSTANT_XRAY
                * z_or_xray
                * cfg.ERG_TO_GEV
                * (cfg.E0 ** 2)
                * area_energy_factor
            )
        else:
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
        prob_source = self.prob_doublet_from_source(z_or_xray, deS)
        prob_dist = self.prob_dist_from_source(
            ra1, de1, sigma1, ra2, de2, sigma2, raS, deS
        )
        return prob_source * prob_dist

    def prob_doublet_background(
        self, de1: float, energy1: float, de2: float, energy2: float
    ) -> float:
        effa1 = self.select_effective_area(de1, energy1) / self.effa_integral
        effa2 = self.select_effective_area(de2, energy2) / self.effa_integral
        return np.cos(de1) * np.cos(de2) * effa1 * effa2 / (4 * np.pi ** 2)

    def set_weights_catalog(self, catalog: catalogs.Catalog) -> None:
        for index in range(len(catalog.names_catalog)):
            if self.flux:
                prob = self.prob_doublet_from_source(
                    catalog.xray_catalog[index], catalog.decs_catalog[index]
                )
            else:
                prob = self.prob_doublet_from_source(
                    catalog.redshifts_catalog[index], catalog.decs_catalog[index]
                )
            self.weights_per_source = np.append(self.weights_per_source, prob)
        self.weights_per_source = self.weights_per_source / np.sum(
            self.weights_per_source
        )

    def select_sources_randomly(
        self, random_generator: np.random.Generator, size: int
    ) -> np.ndarray:
        source_indexes = np.where(
            multinomial.rvs(
                1, self.weights_per_source, random_state=random_generator, size=size
            )
            == 1
        )[1]
        return source_indexes

    def gen_nu_from_source(
        self, z_or_xray: float, dec: float, random_state: np.random.Generator
    ) -> int:
        mu = self.expected_nu_from_source(z_or_xray, dec)
        return poisson.rvs(mu=mu, random_state=random_state)

    def flux_contribute(self, z_or_xray: float, dec: float) -> float:
        if self.flux:
            """
            Given the redshift (or the xray flux) and the declination of a source, determines the contribution
            to the test statistic related to the neutrino flux of the source
            """
            mu = self.expected_nu_from_source(z_or_xray, dec)
            contribute = (
                np.log(0.5) + 2 * np.log(mu) - mu
            )  # Here we assume the limit of low fluxes as valid

        else:
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
        if 90 >= dec > 30:
            effa = self.effective_area_array[cfg.EFFECTIVE_AREA_30_90_DEG_INDEX - 1]
        elif dec <= 30 and dec > 0:
            effa = self.effective_area_array[cfg.EFFECTIVE_AREA_0_30_DEG_INDEX - 1]
        elif dec <= 0 and dec > -5:
            effa = self.effective_area_array[cfg.EFFECTIVE_AREA_MIN5_0_DEG_INDEX - 1]
        elif dec <= -5 and dec > -30:
            effa = self.effective_area_array[
                cfg.EFFECTIVE_AREA_MIN30_MIN5_DEG_INDEX - 1
            ]
        elif dec <= -30 and dec >= -90:
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
        effa1 = self.select_effective_area(dec1, energy1)
        effa2 = self.select_effective_area(dec2, energy2)
        return -np.log(np.cos(dec1) * np.cos(dec2) * effa1 * effa2)

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
        c1 = self.flux_contribute(z_or_xray, decS)
        c2 = self.unc_contribute(sigma1, sigma2)
        c3 = self.dir_contribute(ra1, dec1, ra2, dec2, raS, decS, sigma1, sigma2)
        c4 = self.noise_contribute(dec1, dec2, energy1, energy2)
        ts = c1 + c2 + c3 + c4
        return ts
