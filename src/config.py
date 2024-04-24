ALLOWED_CATALOGS = ["turin", "milliquas"]
ALLOWED_RECONSTRUCTIONS = ["splinempe", "millipede"]
TURIN_INDEX = 0
MILLIQUAS_INDEX = 1
SPLINEMPE_INDEX = 0
MILLIPEDE_INDEX = 1
DEFAULT_CATALOG = ALLOWED_CATALOGS[TURIN_INDEX]
DEFAULT_RECO = ALLOWED_RECONSTRUCTIONS[SPLINEMPE_INDEX]
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
EFFECTIVE_AREA_FILENAME = "Effa_all_streams_gold_bronze.txt"
EFFECTIVE_AREA_ENERGY_BINS_INDEX = 0
EFFECTIVE_AREA_30_90_DEG_INDEX = 1
EFFECTIVE_AREA_0_30_DEG_INDEX = 2
EFFECTIVE_AREA_MIN5_0_DEG_INDEX = 3
EFFECTIVE_AREA_MIN30_MIN5_DEG_INDEX = 4
EFFECTIVE_AREA_MIN90_MIN30_DEG_INDEX = 5
HUBBLE_CONSTANT = 73  # km s^-1 Mpc^-1
MPC_TO_KM = 3.086e19
SPEED_OF_LIGHT = 299792458  # m s^-1
MJD_04102023 = 60132
MJD_GOLDBRONZE_START = 58635
DAYS_IN_YEAR = 365
SECONDS_IN_YEAR = DAYS_IN_YEAR * 24 * 60 * 60
FLUX_NU = 1e36  # s-1 GeV-1 at 100 GeV
CONSTANT_XRAY = 1e-4  # GeV-2
E0 = 100  # GeV
SPLINEMPE_FILENAME = "gcn_notices_gold_bronze.txt"
# IC190819A: Updated GCN Notice never reported
# IC191231A: in GCN Circular reported the day after
# IC200227A: Updated GCN Notice never reported
# IC210503A: Never reconstructed with SplineMPE because detector
# was in a test run configuration. Missing in SplineMPE catalog
# IC210510A: Updated GCN Notice never reported
# IC211117A: Updated GCN Notice never reported
# IC221223A: Updated GCN Notice never reported
# IC230120A: Updated GCN Notice never reported
SPLINEMPE_EXCEPTIONS = [
    "IC191231A",
    "IC190819A",
    "IC200227A",
    "IC210510A",
    "IC211117A",
    "IC221223A",
    "IC230220A",
]
# IC200120A: likely background
# IC230823A: likely background
SPLINEMPE_BACKGROUND = ["IC200120A", "IC230823A"]
SPLINEMPE_GCN_START = "<tr align=left>\n"
SPLINEMPE_INDEX_START = 65
SPLINEMPE_COMMENT_START = "<!--\n"
MILLIPEDE_FILENAME = "IC_Alerts_Table_Energies.csv"
MILLIPEDE_IC_NAME = "IC_NAME"
MILLIPEDE_RA = "RA"
MILLIPEDE_RA_PLUS = "RA_ERR_P"
MILLIPEDE_RA_MINUS = "RA_ERR_M"
MILLIPEDE_DEC = "DEC"
MILLIPEDE_DEC_PLUS = "DEC_ERR_P"
MILLIPEDE_DEC_MINUS = "DEC_ERR_M"
MILLIPEDE_ENERGY = "E [GeV]"
RATIO_90_TO_SIGMA = 2.146
RATIO_68_TO_SIGMA = 1.515
RATIO_50_TO_SIGMA = 1.177
TOTAL_SCRAMBLINGS_SPLINEMPE_TURIN = 300000
TOTAL_SCRAMBLINGS_SPLINEMPE_MILLIQUAS = 70000
TOTAL_SCRAMBLINGS_MILLIPEDE_TURIN = 10000
TOTAL_SCRAMBLINGS_MILLIPEDE_TURIN_XRAY = 100000
TOTAL_SCRAMBLINGS_MILLIPEDE_MILLIQUAS = 500
ROUND_ANGLE = 360  # deg
SPLINEMPE_ANG_DIST_FAST_SELECTION = 4  # deg
MILLIPEDE_ANG_DIST_FAST_SELECTION = 5  # deg
SPLINEMPE_SEARCH_RADIUS = 3  # deg
MILLIPEDE_SEARCH_RADIUS = 4  # deg
TEST_STATISTIC_EMPTY_SCRAMBLE = -1000
ERG_TO_GEV = 624.151
TEV_TO_GEV = 1e3
TEST_STATISTIC_FILENAME = "test_statistic"
FLUX_CHOICES = ["False", "True"]
TRUE_INDEX = 1
FALSE_INDEX = 0

# Flux plot
LINEWIDTH_FLUXPLOT = 3
FLUX1 = 1e36
FLUX2 = 1e39
FLUX3 = 1e34
FLUX4 = 1e42
DECS = 8.5

# plot sed
E_NU_2022 = 1.8399e05  # GeV
E_NU_2023 = 1.2729e05  # GeV
GEV_TO_ERG = 1.60218e-3
TEV_TO_ERG = 1e3 * GEV_TO_ERG
NUMBER_OF_NU = 2
SUPERIOR_ENERGY_BOUND = 2e7  # GeV
M2_TO_CM2 = 1e4
FLUX_NU_MIN_1068 = 2.9e-11 * 1.6  # erg*cm-2*s-1, at 1 TeV
FLUX_NU_MAX_1068 = 7.1e-11 * 1.6  # erg*cm-2*s-1, at 1 TeV
GAMMA_MIN_1068 = 3.0  # spectral index
GAMMA_MAX_1068 = 3.4  # spectral index
MIN_E_NU_1068 = 1.5e3  # GeV
MAX_E_NU_1068 = 1.5e4  # GeV
MAX_E_7YRS = 1e6  # TeV
MIN_E_7YRS = 1e-1  # TeV
LOG_WIDTH_7YRS_BINS = 0.5
INTRINSIC_HARD_XRAY_FLUX_1068 = 206.00e-12  # erg*cm-2*s-1
INTRINSIC_HARD_XRAY_FLUX_ERR_1068 = 22e-12
INTRINSIC_GAMMA_1068 = 2.37
INTRINSIC_GAMMA_ERR_1068 = 0.1
INTRINSIC_HARD_XRAY_FLUX_7469 = 64.30e-12  # erg*cm-2*s-1
INTRINSIC_HARD_XRAY_FLUX_ERR_7469 = 8e-12
INTRINSIC_GAMMA_7469 = 2.09
INTRINSIC_GAMMA_ERR_7469 = 0.12
ENERGY1_INTRINSIC_FLUX = 30  # keV
ENERGY2_INTRINSIC_FLUX = 70  # keV

# coincidence plots
RA_INDEX = 0
DE_INDEX = 1
PLUS_INDEX = 1
MINUS_INDEX = 0
CIRCULAR_23_COORDS = [345.85, 9.14]
CIRCULAR_23_ERRS = [[1.04, 0.90], [0.76, 0.81]]
CIRCULAR_22_COORDS = [346.11, 8.91]
CIRCULAR_22_ERRS = [[1.33, 1.26], [1.01, 0.95]]
NGC7469_COORDS = [345.8156, 8.873861]
CIRCULAR_PLOT_FILENAME = "GCNCirculars_ngc7469"
NOTICE_22_COORDS = [345.757, 8.856]
NOTICE_22_ERR = 0.66
NOTICE_23_COORDS = [345.825, 9.007]
NOTICE_23_ERR = 0.513
NOTICE_PLOT_FILENAME = "GCNNotices_ngc7469"

# plot test statistic
SPLINEMPE = "splinempe"
MILLIPEDE = "millipede"
TURIN = "turin"
MILLIQUAS = "milliquas"
REDSHIFT = "redshift"
XRAY = "xray"
TS_FILENAME = "test_statistic"
FIGSIZE_TS = (5, 3.5)
NBINS_TS = 50
HISTTYPE_TS = "step"
LINEWIDTH_TS = 2
TURIN_LABEL_TS = "Background distrubution,\nTurin catalog"
MILLIQUAS_LABEL_TS = "Background distrubution,\nMilliquas catalog"
AXVCOLOR_TS = "red"
AXVLINESTYLE_TS = "--"
AXVLABEL_TS = "IC220424A & IC230416A\nwith NGC 7469"
XLABEL_TS = "Negative test statistic"
YLABEL_TS = "N Datasets"
TITLE_TS = "Goodness of Fit test"
XRAYTITLE_TS = "Goodness of Fit test, X-ray weighting"
MILLIPEDETITLE_TS = "Goodness of Fit test, Millipede"
MILLIPEDEXRAYTITLE_TS = "Goodness of Fit test, Millipede, X-ray weighting"
FONTSIZE_TS = "large"
AXISSCALE_TS = "log"
FIGNAME_TS = "TS_distr"
BBOX_INCHES = "tight"