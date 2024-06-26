SPLINEMPE = "splinempe"
MILLIPEDE = "millipede"

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
TOTAL_SCRAMBLINGS_SPLINEMPE_INDEX: int = 0
TOTAL_SCRAMBLINGS_MILLIPEDE_INDEX: int = 1
ROUND_ANGLE = 360  # deg
SPLINEMPE_ANG_DIST_FAST_SELECTION = 6.  # deg
MILLIPEDE_ANG_DIST_FAST_SELECTION = 5.  # deg
SPLINEMPE_SEARCH_RADIUS = 5.  # deg
MILLIPEDE_SEARCH_RADIUS = 4.  # deg
TEV_TO_GEV = 1e3