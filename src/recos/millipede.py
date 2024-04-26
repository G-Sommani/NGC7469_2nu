from typing import Final
import pandas as pd  # type: ignore
import numpy as np
from pathlib import Path

from . import Reco
from . import config_recos as cfg


class Millipede(Reco):
    def __init__(self) -> None:
        super().__init__()
        self.RAs_ERR_PLUS: np.ndarray = np.array([])
        self.RAs_ERR_MINUS: np.ndarray = np.array([])
        self.DECs_ERR_PLUS: np.ndarray = np.array([])
        self.DECs_ERR_MINUS: np.ndarray = np.array([])
        self.ang_dist_fast_selection = cfg.MILLIPEDE_ANG_DIST_FAST_SELECTION
        self.search_radius = cfg.MILLIPEDE_SEARCH_RADIUS
        self.total_scramblings_index = cfg.TOTAL_SCRAMBLINGS_MILLIPEDE_INDEX
        self.reco_name = cfg.MILLIPEDE

    def millipede_area(self, index: int) -> float:
        """
        Given the index of the alert, returns the angular area
        """
        phi1_deg = self.RAs[index] - self.RAs_ERR_MINUS[index]
        phi2_deg = self.RAs[index] + self.RAs_ERR_PLUS[index]
        delta1_deg = self.DECs[index] - self.DECs_ERR_MINUS[index]
        delta2_deg = self.DECs[index] + self.DECs_ERR_PLUS[index]
        phi1_rad = np.deg2rad(phi1_deg)
        phi2_rad = np.deg2rad(phi2_deg)
        delta1_rad = np.deg2rad(delta1_deg)
        delta2_rad = np.deg2rad(delta2_deg)
        A = (phi2_rad - phi1_rad) * (np.sin(delta2_rad) - np.sin(delta1_rad))
        return A

    def load_reco_data(self, data_path: Path) -> None:
        alerts_df = pd.read_csv(data_path / cfg.MILLIPEDE_FILENAME)
        self.RAs = alerts_df[cfg.MILLIPEDE_RA].to_numpy()
        self.DECs = alerts_df[cfg.MILLIPEDE_DEC].to_numpy()
        self.RAs_ERR_PLUS = alerts_df[cfg.MILLIPEDE_RA_PLUS].to_numpy()
        self.DECs_ERR_PLUS = alerts_df[cfg.MILLIPEDE_DEC_PLUS].to_numpy()
        self.RAs_ERR_MINUS = alerts_df[cfg.MILLIPEDE_RA_MINUS].to_numpy()
        self.DECs_ERR_MINUS = alerts_df[cfg.MILLIPEDE_DEC_MINUS].to_numpy()
        self.NAMEs = alerts_df[cfg.MILLIPEDE_IC_NAME].to_numpy()
        self.ENERGIES = alerts_df[cfg.MILLIPEDE_ENERGY].to_numpy()
        for i in range(len(alerts_df)):
            area = self.millipede_area(i)
            sigma = np.sqrt(area / np.pi) / cfg.RATIO_90_TO_SIGMA
            self.sigmas = np.append(self.sigmas, sigma)


RECO_CLASS: Final[type[Reco]] = Millipede
