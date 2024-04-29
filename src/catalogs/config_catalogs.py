# catalogs
TURIN = "turin"
MILLIQUAS = "milliquas"

MILLIQUAS_FILENAME = "milliquas.txt"
MILLIQUAS_URL = "https://quasars.org/milliquas.zip"
MILLIQUAS_ZIP = "milliquas.zip"
MILLIQUAS_RA = "RA"
MILLIQUAS_DEC = "DEC"
MILLIQUAS_Z = "Z"
MILLIQUAS_TYPE = "Type"
MILLIQUAS_NAME = "Name"
MILLIQUAS_COLSPECS = [
    (0, 11),
    (12, 23),
    (25, 50),
    (51, 55),
    (56, 61),
    (62, 67),
    (68, 70),
    (71, 72),
    (73, 74),
    (76, 82),
    (83, 88),
    (90, 95),
    (97, 99),
    (101, 103),
    (105, 126),
    (128, 149),
    (151, 172),
    (174, 195),
]
MILLIQUAS_HEADER = [
    MILLIQUAS_RA,
    MILLIQUAS_DEC,
    MILLIQUAS_NAME,
    MILLIQUAS_TYPE,
    "Rmag",
    "Bmag",
    "Comment",
    "R",
    "B",
    MILLIQUAS_Z,
    "Cite",
    "Zcite",
    "RXpct",
    "Qpct",
    "Xname",
    "Rname",
    "Lobe1",
    "Lobe2",
]
MILLIQUAS_Z_CUT = 0.5
MILLIQUAS_AGN_CATEGORIES = [
    "ARX",
    "AX",
    "QRX",
    "QX",
    "BRX",
    "BX",
    "KRX",
    "KX",
    "NRX",
    "NX",
]
TURIN_FILENAME = "Turin_Catalogue_Table2_SourceProperties.txt"
TURIN_NAMES_FILENAME = "Turin_Catalogue_Table1_SourceNames.txt"
BASS_XRAY_FILENAME = "Bass_Catalogue_SourceXRay.txt"
TURIN_URL = "https://cdsarc.cds.unistra.fr/ftp/J/A+A/659/A32/tablea2.dat"
TURIN_NAMES_URL = "https://cdsarc.cds.unistra.fr/ftp/J/A+A/659/A32/tablea1.dat"
BASS_XRAY_URL = "https://content.cld.iop.org/journals/0067-0049/233/2/17/revision1/apjsaa96adt12_mrt.txt"
TURIN_ID = "ID"
TURIN_SYCAT = "SYCAT"
TURIN_SWIFT = "BAT105-Swift"
BASS_SWIFT = "SwiftID"
XRAY_FLUX_NAME = "F14-195-intr"
TURIN_WISE = "WISE"
TURIN_NAME_SOURCE = "IBIS4CAT-IGR"
TURIN_RAh = "RAh"
TURIN_RAm = "RAm"
TURIN_RAs = "RAs"
TURIN_DEd = "DEd"
TURIN_DEm = "DEm"
TURIN_DEs = "DEs"
TURIN_z = "z"
BASS_SKIPROWS = 28
BASS_SEP = "\s+"
FLUX_FACTOR = 1e-12
TURIN_HEADER = [
    TURIN_ID,
    "SYCAT",
    "WISE",
    TURIN_RAh,
    TURIN_RAm,
    TURIN_RAs,
    TURIN_DEd,
    TURIN_DEm,
    TURIN_DEs,
    TURIN_z,
    "LD",
    "Class",
    "r_z",
    "ClassL",
    "r_ClassL1",
    "r_ClassL2",
]
TURIN_NAMES_HEADER = [
    TURIN_ID,
    TURIN_SYCAT,
    "SUMSS",
    "NVSS",
    TURIN_WISE,
    "Pan-STARRS",
    "ROSAT",
    "3PBC",
    TURIN_SWIFT,
    TURIN_NAME_SOURCE,
]
BASS_XRAY_HEADER = [
    BASS_SWIFT,
    "F2-10-obs",
    "F14-195-obs",
    "F2-10-intr",
    "F20-50-intr",
    "F14-150-intr",
    XRAY_FLUX_NAME,
]
TURIN_NAMES_COLSPECS = [
    (0, 2),
    (4, 18),
    (19, 38),
    (40, 58),
    (60, 79),
    (80, 105),
    (107, 122),
    (124, 140),
    (142, 161),
    (162, 182),
]
TURIN_HEADER_PRESENT = None
TOTAL_SCRAMBLINGS_SPLINEMPE_TURIN = 300000
TOTAL_SCRAMBLINGS_SPLINEMPE_TURIN_XRAY = 300000
TOTAL_SCRAMBLINGS_SPLINEMPE_MILLIQUAS = 70000
TOTAL_SCRAMBLINGS_MILLIPEDE_TURIN = 10000
TOTAL_SCRAMBLINGS_MILLIPEDE_TURIN_XRAY = 100000
TOTAL_SCRAMBLINGS_MILLIPEDE_MILLIQUAS = 500