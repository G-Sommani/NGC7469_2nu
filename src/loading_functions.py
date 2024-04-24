import os
from pathlib import Path
from typing import Tuple


def define_paths() -> Tuple[Path, Path]:

    print("Definition of paths...")

    # Definition of paths
    cwd = Path(os.getcwd())
    data_path = cwd / "../data"
    data_results_path = cwd / "../data_results"

    return data_path, data_results_path
