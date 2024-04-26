from typing import Final

from . import Reco


class SplineMPE(Reco):
    def __init__(self) -> None:
        super().__init__()


RECO_CLASS: Final[type[Reco]] = SplineMPE
