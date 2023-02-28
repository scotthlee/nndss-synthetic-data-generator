"""
A collection of NETSS-related constants.

NETSS: National Electronic Telecommunications System for Surveillance
"""


# name strings for the variables in the synthetic NETSS data
NAME_AGE      = 'AGE'
NAME_SEX      = 'SEX'
NAME_RACE     = 'RACE'
NAME_HISPANIC = 'HISPANIC'
NAME_COUNTY   = 'COUNTY'
NAME_CASSTAT  = 'CASSTAT'

# all names in a set
NAME_STRING_SET = {
    NAME_AGE,
    NAME_SEX,
    NAME_RACE,
    NAME_HISPANIC,
    NAME_COUNTY,
    NAME_CASSTAT
}

# age: remap 999 (unknown) to -1, to span contiguous range [-1, 120];
# fractions will be converted to 0
AGE_UNKNOWN       = 999
AGE_UNKNOWN_REMAP = -1
AGE_MAP = {
    AGE_UNKNOWN:AGE_UNKNOWN_REMAP
}

AGETYPE_YEARS     = 0
AGETYPE_UNKNOWN   = 9

# sex: remap 9 (unknown) to 0, to span contiguous range [0, 2]
#          M    F    U
SEX_MALE          = 1
SEX_FEMALE        = 2
SEX_UNKNOWN       = 9
SEX_UNKNOWN_REMAP = 0
SEX_MAP = {
    SEX_MALE   :1,
    SEX_FEMALE :2,
    SEX_UNKNOWN:SEX_UNKNOWN_REMAP
}

# race: remap as above to span contiguous range [0, 5]
#               NA,  A,   B,   W,   O,   U
RACE_NATIVE_AMERICAN = 1
RACE_ASIAN           = 2
RACE_BLACK           = 3
RACE_WHITE           = 5
RACE_OTHER           = 8
RACE_UNKNOWN         = 9
RACE_UNKNOWN_REMAP   = 0
RACE_MAP = {
    RACE_NATIVE_AMERICAN:1,
    RACE_ASIAN:2,
    RACE_BLACK:3,
    RACE_WHITE:4,
    RACE_OTHER:5,
    RACE_UNKNOWN:RACE_UNKNOWN_REMAP
}

# hispanic: remap 9 (unknown) to 0, to span contiguous range [0, 2]
#               Y    N    U
HISPANIC_UNKNOWN       = 9
HISPANIC_UNKNOWN_REMAP = 0
HISPANIC_YES           = 1
HISPANIC_NO            = 2
HISPANIC_MAP = {
    HISPANIC_YES    :1,
    HISPANIC_NO     :2, 
    HISPANIC_UNKNOWN:HISPANIC_UNKNOWN_REMAP
}

# county: this map and its inverse will be built below
COUNTY_UNKNOWN       = 999
COUNTY_UNKNOWN_REMAP = 0

# casstat: remap 9 (unknown) to 0
#              C    P    S    U
CASSTAT_UNKNOWN       = 9
CASSTAT_UNKNOWN_REMAP = 0
CASSTAT_CONFIRMED     = 1
CASSTAT_PROBABLE      = 2
CASSTAT_SUSPECT       = 3
CASSTAT_MAP = {
    CASSTAT_CONFIRMED:1,
    CASSTAT_PROBABLE :2,
    CASSTAT_SUSPECT  :3,
    CASSTAT_UNKNOWN:CASSTAT_UNKNOWN_REMAP
}

# min and max bin values for each variable;
# values are at the LEFT EDGE of the bin
AGE_MIN      = -1
AGE_MAX      = 120
SEX_MIN      = 0
SEX_MAX      = 2
RACE_MIN     = 0
RACE_MAX     = 5
HISPANIC_MIN = 0
HISPANIC_MAX = 2
CASSTAT_MIN  = 0
CASSTAT_MAX  = 3
COUNTY_MIN   = 0
# county_max depends on the data

INV_AGE_MAP      = {v:k for k,v in AGE_MAP.items()}
INV_SEX_MAP      = {v:k for k,v in SEX_MAP.items()}
INV_RACE_MAP     = {v:k for k,v in RACE_MAP.items()}
INV_HISPANIC_MAP = {v:k for k,v in HISPANIC_MAP.items()}
INV_CASSTAT_MAP  = {v:k for k,v in CASSTAT_MAP.items()}
