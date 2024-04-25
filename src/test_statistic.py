from loading_functions import Loader
import config as cfg
import numpy as np


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

    def select_effective_area(self, dec: float, energy: float):
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
