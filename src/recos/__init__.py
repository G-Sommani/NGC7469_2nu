from abc import ABC
import importlib


class Reco(ABC):
    def __init__(self) -> None:
        pass


def initiate_reco(reco_name: str) -> Reco:
    module = importlib.import_module(f"{__name__}.{reco_name.lower()}")
    return module.RECO_CLASS()
