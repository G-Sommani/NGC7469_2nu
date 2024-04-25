from abc import ABC
import importlib


class Catalog(ABC):

    catalog_name: str = __name__
    xray: bool = False
    zipname_data: str | None = None
    filename_names: str | None = None
    url_names: str | None = None
    filename_xray: str | None = None
    url_xray: str | None = None
    filename_data: str
    url_data: str

    def __init__(self, xray: bool = False) -> None:
        self.xray = xray
        pass


def initiate_catalog(catalog_name: str, xray: bool = False) -> Catalog:
    module = importlib.import_module(f"{__name__}.{catalog_name.lower()}")
    return module.CATALOG_CLASS(xray=xray)
