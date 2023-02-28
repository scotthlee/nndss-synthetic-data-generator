"""
A collection of HL7-related constants.
"""


# name strings for the variables in the synthetic HL7 data
NAME_AGE         = 'AGE'
NAME_SEX         = 'SEX'
NAME_RACE        = 'RACE'
NAME_ETHNICITY   = 'ETHNICITY'
NAME_COUNTY      = 'COUNTY'
NAME_CASE_STATUS = 'CASE_STATUS'
NAME_PREGNANT    = 'PREGNANT'

# all names in a set
NAME_STRING_SET = {
    NAME_AGE,
    NAME_SEX,
    NAME_RACE,
    NAME_ETHNICITY,
    NAME_COUNTY,
    NAME_CASE_STATUS,
    NAME_PREGNANT
}

# age: remap 9999 (unknown) to -1, to span contiguous range [-1, 120];
# fractions will be converted to 0
AGE_UNKNOWN       = 9999
AGE_UNKNOWN_REMAP = -1
AGE_MAP = {
    AGE_UNKNOWN:AGE_UNKNOWN_REMAP
}

# unknown date value in some HL7 fields
DATE_UNKNOWN      = '99999999'

AGE_UNITS_YEARS   = 'a'
AGE_UNITS_MONTHS  = 'mo'
AGE_UNITS_WEEKS   = 'wk'
AGE_UNITS_DAYS    = 'd'
AGE_UNITS_OTHER   = 'oth' 
AGE_UNITS_UNKNOWN = 'unk'

# age units for which an age can be computed
COMPUTABLE_AGE_UNITS = {
    AGE_UNITS_YEARS, AGE_UNITS_MONTHS, AGE_UNITS_WEEKS, AGE_UNITS_DAYS
}

# sex:
SEX_UNKNOWN = 'u'
SEX_MALE    = 'm'
SEX_FEMALE  = 'f'
SEX_VALUE_UNKNOWN = 0.0
SEX_VALUE_MALE    = 1.0
SEX_VALUE_FEMALE  = 2.0
SEX_MAP = {
    SEX_UNKNOWN:SEX_VALUE_UNKNOWN,
    SEX_MALE:SEX_VALUE_MALE,
    SEX_FEMALE:SEX_VALUE_FEMALE
}

SEX_MAP_VALUES = {v for v in SEX_MAP.values()}

# race
RACE_UNKNOWN           = 'nullflavor'
RACE_UNKNOWN_REMAP     = 0
RACE_UNKNOWN_REMAP_STR = '0'

# ethnicity
ETHNICITY_UNKNOWN      = 'unknown'
ETHNICITY_HISPANIC     = 'hispanic or latino'
ETHNICITY_NON_HISPANIC = 'not hispanic or latino'
ETHNICITY_OTHER        = 'other'
ETHNICITY_MAP = {
    ETHNICITY_UNKNOWN:0.0,
    ETHNICITY_HISPANIC:1.0,
    ETHNICITY_NON_HISPANIC:2.0,
    ETHNICITY_OTHER:3.0
}

ETHNICITY_MAP_VALUES = {v for v in ETHNICITY_MAP.values()}

# case-status
CASE_STATUS_UNKNOWN   = 'unknown'
CASE_STATUS_NOT       = 'not a case'
CASE_STATUS_SUSPECTED = 'suspected'
CASE_STATUS_PROBABLE  = 'probable diagnosis'
CASE_STATUS_CONFIRMED = 'confirmed present'
CASE_STATUS_VALUE_UNKNOWN   = 0.0
CASE_STATUS_VALUE_NOT       = 1.0
CASE_STATUS_VALUE_SUSPECTED = 2.0
CASE_STATUS_VALUE_PROBABLE  = 3.0
CASE_STATUS_VALUE_CONFIRMED = 4.0
CASE_STATUS_MAP = {
    CASE_STATUS_UNKNOWN:CASE_STATUS_VALUE_UNKNOWN,
    CASE_STATUS_NOT:CASE_STATUS_VALUE_NOT,
    CASE_STATUS_SUSPECTED:CASE_STATUS_VALUE_SUSPECTED,
    CASE_STATUS_PROBABLE:CASE_STATUS_VALUE_PROBABLE,
    CASE_STATUS_CONFIRMED:CASE_STATUS_VALUE_CONFIRMED
}

CASE_STATUS_MAP_VALUES = {v for v in CASE_STATUS_MAP.values()}

# county: this map and its inverse will be built below
COUNTY_UNKNOWN           = 'unk'
COUNTY_UNKNOWN_REMAP     = 0
COUNTY_UNKNOWN_REMAP_STR = '0'

# pregnant
PREGNANT_UNKNOWN = 'unk'
PREGNANT_NO      = 'n'
PREGNANT_YES     = 'y'
PREGNANT_VALUE_UNKNOWN = 0.0
PREGNANT_VALUE_NO      = 1.0
PREGNANT_VALUE_YES     = 2.0
PREGNANT_MAP = {
    PREGNANT_UNKNOWN:PREGNANT_VALUE_UNKNOWN,
    PREGNANT_NO:PREGNANT_VALUE_NO,
    PREGNANT_YES:PREGNANT_VALUE_YES
}

PREGNANT_MAP_VALUES = {v for v in PREGNANT_MAP.values()}

# min and max bin values for each variable;
# values are at the LEFT EDGE of the bin
AGE_MIN        = -1
AGE_MIN_FLOAT  = -1.0
AGE_MAX        = 120
AGE_MAX_FLOAT  = 120.0
SEX_MIN        = 0
SEX_MIN_FLOAT  = 0.0
SEX_MAX        = 2
SEX_MAX_FLOAT  = 2.0
RACE_MIN       = 0
RACE_MIN_FLOAT = 0.0
ETHNICITY_MIN       = 0
ETHNICITY_MIN_FLOAT = 0.0
ETHNICITY_MAX       = 3
ETHNICITY_MAX_FLOAT = 3.0
CASE_STATUS_MIN        = 0
CASE_STATUS_MIN_FLOAT  = 0.0
CASE_STATUS_MAX        = 4
CASE_STATUS_MAX_FLOAT  = 4.0
COUNTY_MIN       = 0
COUNTY_MIN_FLOAT = 0.0
PREGNANT_MIN       = 0
PREGNANT_MIN_FLOAT = 0.0
PREGNANT_MAX       = 2
PREGNANT_MAX_FLOAT = 2.0

# race_max depends on the data
# county_max depends on the data

# inverse maps
INV_AGE_MAP = {v:k for k,v in AGE_MAP.items()}
INV_SEX_MAP = {v:k for k,v in SEX_MAP.items()}
INV_ETHNICITY_MAP = {v:k for k,v in ETHNICITY_MAP.items()}
INV_CASE_STATUS_MAP = {v:k for k,v in CASE_STATUS_MAP.items()}
INV_PREGNANT_MAP = {v:k for k,v in PREGNANT_MAP.items()}

# notif_result_status
NOTIF_RESULT_STATUS_FINAL        = 'f'
NOTIF_RESULT_STATUS_CORRECTION   = 'c'
NOTIF_RESULT_STATUS_NOT_OBTAINED = 'x'

# Fields for the output record, which are either the JSON field name strings
# or the CSV header strings.
# These must be in the same order as the preprocessed HL7 data fields.
FIELD_REPORT_DATE = 'report_dt_new'          # 0
FIELD_COUNT       = 'count'                  # 1
FIELD_AGE         = 'age'                    # 2
FIELD_AGE_UNITS   = 'age_units'              # 3
FIELD_SEX         = 'sex'                    # 4
FIELD_ETHNICITY   = 'ethnicity_txt'          # 5
FIELD_RACE        = 'race_mapped'            # 6
FIELD_CASE_STATUS = 'case_status_txt'        # 7
FIELD_BIRTH_DATE  = 'birth_date_str'         # 8
FIELD_NRS         = 'notif_result_status'    # 9
FIELD_PREGNANT    = 'pregnant'               # 10
FIELD_COUNTY      = 'report_county'          # 11
FIELD_DATE1       = 'first_elec_submit_dt'   # 12
FIELD_SUBJ_COUNTY = 'subj_county'            # 13
FIELD_DATE2       = 'diag_dt'                # 14
FIELD_DATE3       = 'died_dt'                # 15
FIELD_DATE4       = 'earliest_cnty_dt'       # 16
FIELD_DATE5       = 'earliest_state_dt'      # 17
FIELD_DATE6       = 'hosp_admit_dt'          # 18
FIELD_DATE7       = 'illness_onset_dt'       # 19
FIELD_DATE8       = 'invest_start_dt'        # 20
FIELD_DATE9       = 'phd_notif_dt'           # 21

# Update these offsets and the empty line template if file header changes.
OFFSET_REPORT_DT_NEW        =  0
OFFSET_COUNT                =  1
OFFSET_FIRST_ELEC_SUBMIT_DT = 12

# The date used by the preprocessor to sort the data. Also used as the
# date tuple anchor date.
FIELD_REF_DATE    = FIELD_DATE1
OFFSET_REF_DATE   = OFFSET_FIRST_ELEC_SUBMIT_DT

# empty line template, used in _merge_files.
EMPTY_LINE = ',0,,,,,,,,,,,{0},,,,,,,,,\n'

DATE_TUPLE_FIELDS = [
    FIELD_REPORT_DATE, # 'report_dt_new'
    FIELD_DATE1,       # 'first_elec_submit_dt'
    FIELD_DATE2,       # 'diag_dt'
    FIELD_DATE3,       # 'died_dt'
    FIELD_DATE6,       # 'hosp_admit_dt',
    FIELD_DATE7,       # 'illness_onset_dt',
    FIELD_DATE8,       # 'invest_start_dt'
]

OUTPUT_FIELDS = [
    FIELD_REPORT_DATE,
    FIELD_COUNT,
    FIELD_AGE,
    FIELD_AGE_UNITS,
    FIELD_SEX,
    FIELD_ETHNICITY,
    FIELD_RACE,
    FIELD_CASE_STATUS,
    FIELD_BIRTH_DATE,
    FIELD_NRS,
    FIELD_PREGNANT,
    FIELD_COUNTY,
    FIELD_DATE1,
    FIELD_SUBJ_COUNTY,
    FIELD_DATE2,
    FIELD_DATE3,
    FIELD_DATE4,
    FIELD_DATE5,
    FIELD_DATE6,
    FIELD_DATE7,
    FIELD_DATE8,
    FIELD_DATE9,
]

# pseudoperson categories (first letter is SEX, second is CASE_STATUS)
PSEUDOPERSON_VALUE_UU =  0.0  # unknown, unknown
PSEUDOPERSON_VALUE_US = 32.0  # unknown, suspected
PSEUDOPERSON_VALUE_UP = 48.0  # unknown, probable
PSEUDOPERSON_VALUE_UC = 64.0  # unknown, confirmed
PSEUDOPERSON_VALUE_MU =  1.0  # male,    unknown
PSEUDOPERSON_VALUE_MS = 33.0  # male,    suspected
PSEUDOPERSON_VALUE_MP = 49.0  # male,    probable
PSEUDOPERSON_VALUE_MC = 65.0  # male,    confirmed
PSEUDOPERSON_VALUE_FU =  2.0  # female,  unknown
PSEUDOPERSON_VALUE_FS = 34.0  # female,  suspected
PSEUDOPERSON_VALUE_FP = 50.0  # female,  probable
PSEUDOPERSON_VALUE_FC = 66.0  # female,  confirmed
PSEUDOPERSON_VALUE_O  = 100.0 # other, values in PSEUDOPERSON_SET_OTHER

PSEUDOPERSON_SET_OTHER = {
    PSEUDOPERSON_VALUE_UC, PSEUDOPERSON_VALUE_UP, PSEUDOPERSON_VALUE_US,
    PSEUDOPERSON_VALUE_MS, PSEUDOPERSON_VALUE_FS
}

PSEUDOPERSON_SYMBOL_MAP = {
    PSEUDOPERSON_VALUE_UU:'UU',
    PSEUDOPERSON_VALUE_US:'US',
    PSEUDOPERSON_VALUE_UP:'UP',
    PSEUDOPERSON_VALUE_UC:'UC',
    PSEUDOPERSON_VALUE_MU:'MU',
    PSEUDOPERSON_VALUE_MS:'MS',
    PSEUDOPERSON_VALUE_MP:'MP',
    PSEUDOPERSON_VALUE_MC:'MC',
    PSEUDOPERSON_VALUE_FU:'FU',
    PSEUDOPERSON_VALUE_FS:'FS',
    PSEUDOPERSON_VALUE_FP:'FP',
    PSEUDOPERSON_VALUE_FC:'FC',
    PSEUDOPERSON_VALUE_O:'OTHER',
}
