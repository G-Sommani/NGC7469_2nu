from abc import ABC, abstractmethod
import importlib
import numpy as np
from pathlib import Path


class Reco(ABC):

    RAs: np.ndarray = np.array([])
    DECs: np.ndarray = np.array([])
    sigmas: np.ndarray = np.array([])
    NAMEs: np.ndarray = np.array([])
    ENERGIES: np.ndarray = np.array([])
    ang_dist_fast_selection: float
    search_radius: float
    total_scramblings_index: int
    reco_name: str

    def __init__(self) -> None:
        pass

    @abstractmethod
    def load_reco_data(self, data_path: Path) -> None:
        pass


def initiate_reco(reco_name: str) -> Reco:
    print(f"Retrieving the alerts reconstructed with {reco_name}...")
    module = importlib.import_module(f"{__name__}.{reco_name.lower()}")
    return module.RECO_CLASS()
