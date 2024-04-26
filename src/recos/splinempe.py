from typing import Final
import numpy as np
from pathlib import Path

from . import Reco
from . import config_recos as cfg


class SplineMPE(Reco):
    def __init__(self) -> None:
        super().__init__()

    def load_reco_data(self, data_path: Path) -> None:
        rev0 = False
        rev1 = False
        has_rev1 = False
        is_data = False
        rev1_names = np.array([])
        splinempe_f = open(data_path / cfg.SPLINEMPE_FILENAME)
        notice_line_index = 0
        for index, line in enumerate(splinempe_f):
            if line == cfg.SPLINEMPE_GCN_START and index > cfg.SPLINEMPE_INDEX_START:
                notice_line_index = 0
                if index < 100:
                    is_data = True
            if line == cfg.SPLINEMPE_COMMENT_START:
                is_data = False
            if is_data:
                if notice_line_index == 2:
                    rev_number = line.split(">")[1].split("<")[0]
                    if rev_number == "1" or rev_number == "2":
                        rev0 = False
                        rev1 = True
                    if rev_number == "0":
                        rev0 = True
                        rev1 = False
                notice_line_index += 1
                if notice_line_index == 4:
                    date = line.split(">")[1].split("<")[0]
                    year = date.split("/")[0]
                    month = date.split("/")[1]
                    day = date.split("/")[2]
                    name = f"IC{year}{month}{day}A"
                    if rev1:
                        if len(rev1_names) > 0:
                            if name in rev1_names:
                                if name == f"IC{year}{month}{day}B":
                                    continue
                                rev1_names[
                                    np.where(rev1_names == name)
                                ] = f"IC{year}{month}{day}B"
                        rev1_names = np.append(rev1_names, name)
                    elif rev0:
                        if len(self.NAMEs) > 0:
                            if name in self.NAMEs:
                                self.NAMEs[
                                    np.where(self.NAMEs == name)
                                ] = f"IC{year}{month}{day}B"
                        if (
                            name in rev1_names or name in cfg.SPLINEMPE_EXCEPTIONS
                        ) and not name in cfg.SPLINEMPE_BACKGROUND:
                            self.NAMEs = np.append(self.NAMEs, name)
                            has_rev1 = True
                        else:
                            has_rev1 = False
                if rev0 and is_data and has_rev1:
                    if notice_line_index == 7:
                        ra = float(line.split(">")[1].split("<")[0])
                        self.RAs = np.append(self.RAs, ra)
                    if notice_line_index == 8:
                        de = float(line.split(">")[1].split("<")[0])
                        self.DECs = np.append(self.DECs, de)
                    if notice_line_index == 10:
                        err_50 = float(line.split(">")[1].split("<")[0]) / 60
                        sigma = np.deg2rad(err_50) / cfg.RATIO_50_TO_SIGMA
                        self.sigmas = np.append(self.sigmas, sigma)
                    if notice_line_index == 11:
                        energy = (
                            float(line.split(">")[1].split("<")[0]) * cfg.TEV_TO_GEV
                        )
                        self.ENERGIES = np.append(self.ENERGIES, energy)


RECO_CLASS: Final[type[Reco]] = SplineMPE
