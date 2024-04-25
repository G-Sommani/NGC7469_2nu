from loading_functions import Loader
import config as cfg
import numpy as np


class TestStatistic:
    def __init__(self) -> None:

        print("Defining the test statistic...")

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
