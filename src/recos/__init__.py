from abc import ABC
import importlib
import numpy as np


class Reco(ABC):

    RAs: np.ndarray = np.array([])
    DECs: np.ndarray = np.array([])
    sigmas: np.ndarray = np.array([])
    NAMEs: np.ndarray = np.array([])
    ENERGIES: np.ndarray = np.array([])

    def __init__(self) -> None:
        pass


def initiate_reco(reco_name: str) -> Reco:
    print(f"Retrieving the alerts reconstructed with {reco_name}...")
    module = importlib.import_module(f"{__name__}.{reco_name.lower()}")
    return module.RECO_CLASS()
