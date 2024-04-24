import os
from pathlib import Path
from typing import List


def define_paths(
    data: bool = True,
    results: bool = False,
    figures: bool = False,
) -> List[Path]:

    print("Definition of paths...")

    # Definition of paths
    cwd = Path(os.getcwd())
    paths = list()
    if data:
        data_path = cwd / "../data"
        paths.append(data_path)
    if results:
        data_results_path = cwd / "../data_results"
        paths.append(data_results_path)
    if figures:
        figures_path = cwd / "../figures"
        paths.append(figures_path)

    return paths
