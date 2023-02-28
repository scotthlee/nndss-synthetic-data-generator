import os
import re
import sys
import copy
import json
import math
import time
import hashlib
import datetime
import numpy as np

from collections import namedtuple
from collections import defaultdict

from scipy import stats

from . import plots_hl7
from . import hl7 as HL7
from . import correlation_matrix

from .ecdf import EmpiricalCDF
from .kernel_density_estimation import KdeECDF

_MISSING_AGE         = 'missing_age_count'
_MISSING_AGE_UNITS   = 'missing_age_units_count'
_MISSING_SEX         = 'missing_sex_count'
_MISSING_RACE        = 'missing_race_count'
_MISSING_ETHNICITY   = 'missing_ethnicity_count'
_MISSING_COUNTY      = 'missing_county_count'
_MISSING_CASE_STATUS = 'missing_case_status_count'
_MISSING_PREGNANT    = 'missing_pregnant_count'
_INVALID_SEX         = 'invalid_sex_count'
_INVALID_RACE        = 'invalid_race_count'
_INVALID_ETHNICITY   = 'invalid_ethnicity_count'
_INVALID_CASE_STATUS = 'invalid_casstat_count'
_INVALID_PREGNANT    = 'invalid_pregnant_count'
_UNHANDLED_AGE_UNITS = 'unhandled_age_units_count'

_NOT_A_CASE = 'not_a_case_count'
_RESULTS_NOT_OBTAINED = 'results_not_obtained_count'

_HL7_RECORD_FIELDS = [
    'event_date', # match point to signal processing code
    'count',
    'age',
    'sex',
    'race',
    'ethnicity',
    'case_status',
    'county',
    'pregnant',
    'first_elec_submit_dt',
    'diag_dt',
    'died_dt',
    'report_dt_new',
    'hosp_admit_dt',
    'illness_onset_dt',
    'invest_start_dt'
]

_HL7_Record = namedtuple('_HL7_Record', _HL7_RECORD_FIELDS)


# namedtuple used for writing the output files;
# fields default to the empty string
_OutputTuple = namedtuple('_OutputTuple', HL7.OUTPUT_FIELDS)
_OutputTuple.__new__.__defaults__ = ('',) * len(_OutputTuple._fields)

# namedtuple used for pseudopersons
_PSEUDOPERSON_FIELDS = ['tuple_index_map', 'tuple_index_dist', 'tuple_ecdf']
_Pseudoperson = namedtuple('_Pseudoperson', _PSEUDOPERSON_FIELDS)

# HL7 data is extremely sparse prior to 2015-01-01; also, some diagnosis dates
# extend several decades into the past. CDC agrees to use 2015 as the cutoff
# date.
_MIN_DATE = datetime.datetime.strptime('2015-01-01', '%Y-%m-%d')

# anchors for tuple computations
_TUPLE_ANCHOR_EARLIEST = 0
_TUPLE_ANCHOR_MEDIAN   = 1
_TUPLE_ANCHOR = _TUPLE_ANCHOR_EARLIEST

# set to True to display debug info
_TRACE = False

_EMPTY_STRING = ''

# a class to hold all of the state variables
class State:

    def reset(self):
        # clear and initialize all state variables
        
        # preprocessed file header
        self.preprocessed_header_strings = None

        # variable used for file input
        self.nonempty_record_count = 0

        # age
        self.age_data = []
        self.age_dict = defaultdict(int)
        self.ecdf_age = None

        # sex
        self.sex_data = []
        self.sex_dict = defaultdict(int)
        self.ecdf_sex = None

        # race
        self.race_data = []
        self.race_dict = defaultdict(int)
        self.race_map = {}
        self.inv_race_map = {}
        self.race_max = None
        self.ecdf_race = None

        # ethnicity
        self.ethnicity_data = []
        self.ethnicity_dict = defaultdict(int)
        self.ecdf_ethnicity = None

        # case status
        self.case_status_data = []
        self.case_status_dict = defaultdict(int)
        self.ecdf_case_status = None

        # county
        self.county_data = []
        self.county_dict = defaultdict(int)
        self.county_map = {}
        self.inv_county_map = {}
        self.county_max = None
        self.ecdf_county = None

        # pregnant
        self.pregnant_data = []
        self.pregnant_dict = defaultdict(int)
        self.ecdf_pregnant = None

        # pseudoperson state management
        # the key to this map is created by _to_pseudoperson_category() below
        self.pseudoperson_map = {}

    def __init__(self):
        self.reset()

# instance of the state class
_state = State()


###############################################################################
def enable_debug():
    global _TRACE
    _TRACE = True


###############################################################################
def get_preprocessed_header():
    return _state.preprocessed_header_strings


###############################################################################
def _day_diff(str_datetime, ref_datetime):
    """
    Compute the difference in days between two dates. The first arg is a
    datetime string and the second is a datetime object.
    """

    delta = None
    if str_datetime is not None and len(str_datetime) > 0:
        # keep only the YYYY-mm-dd portion
        str_date = str_datetime[:10]
        d = datetime.datetime.strptime(str_date, '%Y-%m-%d')
        delta = (d - ref_datetime).days

    return delta


###############################################################################
def _median_datetime(dt_list):
    """
    Compute the median datetime from a list of datetime objects. Do NOT change
    the order of the items in the list.
    """

    new_list = copy.deepcopy(dt_list)
    new_list = sorted(new_list)
    n = len(new_list)

    # if odd length, take the middle element from the sorted array
    if 1 == n % 2:
        # n//2 is the center element with 0-based indexing
        median_dt = new_list[n//2]
    else:
        # even, so compute the datetime between the two middle datetimes
        dt1 = new_list[n//2 - 1]
        dt2 = new_list[n//2]

        day_diff = (dt2 - dt1).days
        # use integer truncation, to avoid partial days
        median_dt = dt1 + datetime.timedelta(days=day_diff//2)

    return median_dt


###############################################################################
def signal_from_anchor_date(file_data):
    """
    Compute the number of cases per day and anchor the count to either the
    earliest or the median date in each record.
    """

    num_nonzero_counts = 0

    first_elec_submit_dt_map = defaultdict(int)
    diag_dt_map              = defaultdict(int)
    died_dt_map              = defaultdict(int)
    report_dt_new_map        = defaultdict(int)
    hosp_admit_dt_map        = defaultdict(int)
    illness_onset_dt_map     = defaultdict(int)
    invest_start_dt_map      = defaultdict(int)

    # maps a datetime to a list of file_data indices, those records sharing
    # the same anchor datetime
    anchor_map = defaultdict(list)

    dates = []
    max_date = None
    for i,r in enumerate(file_data):

        # candidates for the anchor date in this record
        candidates = []

        if r.count > 0:
            num_nonzero_counts += 1
        else:
            continue

        # r.event_date is always present, and is equal to one of these dates below
        #candidates.append(datetime.datetime.strptime(r.event_date, '%Y-%m-%d'))

        # keep only the YYYY-mm-dd portion of the datetime strings
        if r.first_elec_submit_dt is not None and len(r.first_elec_submit_dt) > 0:
            dt = datetime.datetime.strptime(r.first_elec_submit_dt[:10], '%Y-%m-%d')
            first_elec_submit_dt_map[dt] += 1
            candidates.append(dt)
        if r.diag_dt is not None and len(r.diag_dt) > 0:
            dt = datetime.datetime.strptime(r.diag_dt[:10], '%Y-%m-%d')
            diag_dt_map[dt] += 1
            candidates.append(dt)
        if r.died_dt is not None and len(r.died_dt) > 0:
            dt = datetime.datetime.strptime(r.died_dt[:10], '%Y-%m-%d')
            died_dt_map[dt] += 1
            candidates.append(dt)
        if r.report_dt_new is not None and len(r.report_dt_new) > 0:
            dt = datetime.datetime.strptime(r.report_dt_new[:10], '%Y-%m-%d')
            report_dt_new_map[dt] += 1
            candidates.append(dt)
        if r.hosp_admit_dt is not None and len(r.hosp_admit_dt) > 0:
            dt = datetime.datetime.strptime(r.hosp_admit_dt[:10], '%Y-%m-%d')
            hosp_admit_dt_map[dt] += 1
            candidates.append(dt)
        if r.illness_onset_dt is not None and len(r.illness_onset_dt) > 0:
            dt = datetime.datetime.strptime(r.illness_onset_dt[:10], '%Y-%m-%d')
            illness_onset_dt_map[dt] += 1
            candidates.append(dt)
        if r.invest_start_dt is not None and len(r.invest_start_dt) > 0:
            dt = datetime.datetime.strptime(r.invest_start_dt[:10], '%Y-%m-%d')
            invest_start_dt_map[dt] += 1
            candidates.append(dt)

        # at least one of these dates must exist
        assert len(candidates) > 0
        assert None not in anchor_map

        if _TUPLE_ANCHOR_EARLIEST == _TUPLE_ANCHOR:
            anchor_date = min(candidates)
        else:
            anchor_date = _median_datetime(candidates)

            # if datetime.datetime(2017, 9, 11) == anchor_date:
            #     print('No. candidates: {0}'.format(len(candidates)))
            #     print('\t       report_dt_new: {0}'.format(r.report_dt_new[:10]))
            #     print('\tfirst_elec_submit_dt: {0}'.format(r.first_elec_submit_dt[:10]))
            #     print('\t             diag_dt: {0}'.format(r.diag_dt[:10]))
            #     print('\t             died_dt: {0}'.format(r.died_dt[:10]))
            #     print('\t       hosp_admit_dt: {0}'.format(r.hosp_admit_dt[:10]))
            #     print('\t    illness_onset_dt: {0}'.format(r.illness_onset_dt[:10]))
            #     print('\t     invest_start_dt: {0}'.format(r.invest_start_dt[:10]))

            
        anchor_map[anchor_date].append(i)

    print('Found {0} records with nonzero counts'.format(num_nonzero_counts))

    sorted_dates = sorted([k for k in anchor_map.keys()])
    print('Number of sorted anchor dates: {0}'.format(len(sorted_dates)))
    min_anchor_date = sorted_dates[0]
    max_anchor_date = sorted_dates[-1]
    day_count = (max_anchor_date - min_anchor_date).days
    print('Min anchor date: {0}, max anchor date: {1}, day_count: {2}'.
          format(min_anchor_date, max_anchor_date, day_count))

    # generate YYYY-mm-dd date strings for the entire day range
    dates = []
    cur_date = min_anchor_date
    while cur_date <= max_anchor_date:
        str_date = datetime.datetime.strftime(cur_date, '%Y-%m-%d')
        dates.append(str_date)
        cur_date += datetime.timedelta(days=1)

    #print('len(dates): {0}'.format(len(dates)))
    print('Min date from dates list: {0}'.format(dates[0]))
    print('Max date from dates list: {0}'.format(dates[-1]))
    assert len(sorted_dates) <= len(dates)

    # Construct the signal as the number of nonzero cases on each anchor date.
    # Reference as an offset from the min anchor date.
    signal = [0] * len(dates)

    dt_of_max_signal = None
    max_signal = -1
    for d in sorted_dates:
        date_offset = (d - min_anchor_date).days
        # anchor_map: date (datetime obj) => list of indices of cases sharing
        # that anchor date
        if d in anchor_map:
            cases_this_day = len(anchor_map[d])
            signal[date_offset] = cases_this_day
            if cases_this_day > max_signal:
                max_signal = cases_this_day
                dt_of_max_signal = d

    print('Signal information: ')
    print('\t min signal value: {0}'.format(min(signal)))
    print('\t max signal value: {0}'.format(max(signal)))
    print('\tdate of max value: {0}'.format(dt_of_max_signal))

    # convert the signal to an np.array
    signal = np.array(signal)

    return dates, signal, (first_elec_submit_dt_map, diag_dt_map, died_dt_map, \
        report_dt_new_map, hosp_admit_dt_map, illness_onset_dt_map,           \
        invest_start_dt_map)


###############################################################################
def _to_pseudoperson_category(sex_value, cs_value):
    """
    Convert the given SEX and CASE_STATUS values to a pseudoperson category.
    """

    # quick transformation, generates a pseudoperson map key
    pp_category = cs_value * 16.0 + sex_value
    if pp_category in HL7.PSEUDOPERSON_SET_OTHER:
        return HL7.PSEUDOPERSON_VALUE_O
    else:
        return pp_category


###############################################################################
def _print_top_pseudoperson_tuples(pp_category, n=5):
    """
    Print the n most common tuples for each pseudoperson type.
    """

    pp_obj = _state.pseudoperson_map[pp_category]

    # maps tuple index to the tuple
    index_map = pp_obj.tuple_index_map
    # maps a tuple index to the tuple multiplicity
    index_dist = pp_obj.tuple_index_dist

    #inv_dist_map = {v:k for k,v in index_dist.items()}
    inv_dist = [(count,tup_index) for tup_index,count in index_dist.items()]
    inv_dist = sorted(inv_dist, key=lambda x: x[0], reverse=True)

    num_printed = 0
    for count, tup_index in inv_dist:
        tup = index_map[tup_index]
        print('\t{0},\tcount: {1}'.format(tup, count))
        num_printed += 1
        if num_printed >= n:
            break


###############################################################################
def _to_pseudoperson_obj(pp_category, tuple_dist):
    """
    Use the pseudoperson category and tuple distribution to create a
    pseudoperson map entry.
    """

    symbol = HL7.PSEUDOPERSON_SYMBOL_MAP[pp_category]
    print('pp_category: "{0}", unique_tuple count: {1}, total tuples: {2}'.
          format(symbol, len(tuple_dist), sum(tuple_dist.values())))

    # assign tuple indices and build a tuple index dist for the ECDF
    tuple_index = 0
    # maps a tuple index to a tuple
    tuple_index_map = {}
    # maps a tuple index to a count (multiplicity of that tuple)
    tuple_index_dist = defaultdict(int)
    for tup,count in tuple_dist.items():
        tuple_index_map[tuple_index] = tup
        tuple_index_dist[tuple_index] = count
        tuple_index += 1

    #print('\tMax tuple multiplicity: {0}'.format(max(tuple_index_dist.values())))

    unique_tuple_count = len(tuple_index_map)
    assert unique_tuple_count == len(tuple_index_dist)
    tuple_ecdf = EmpiricalCDF(tuple_index_dist)

    # compute and save the tuple ECDF for this pseudoperson
    pp = _Pseudoperson(
        tuple_index_map = tuple_index_map,
        tuple_index_dist = tuple_index_dist,
        tuple_ecdf = tuple_ecdf
    )

    return pp


###############################################################################
def _compute_pseudoperson_distributions(file_data):
    """
    Pseudopersons are unique combinations of SEX-CASE_STATUS values. The valid
    SEX values are {U, M, F} and the valid CASE_STATUS values are {U, S, P, C}.
    The not-a-case value for CASE_STATUS is ignored, per CDC guidance.

    This function partitions the file_data into pseudopersons and computes
    date tuple empirical cumulative distribution functions for each.
    """

    pseudoperson_map = defaultdict(list)

    # partition the data into pseudoperson categories

    # 'file_data' is a list of _HL7_Record objects
    for i, r in enumerate(file_data):
        count = r.count
        # skip 0-count records, these have no data
        if 0 == count:
            continue

        # sex is one of the floating point values in SEX_MAP_VALUES
        # cs is one of the floating point values in CASE_STATUS_VALUES
        sex_value = r.sex
        cs_value  = r.case_status

        key = _to_pseudoperson_category(sex_value, cs_value)
        pseudoperson_map[key].append(r)

    print('\nPseudoperson distribution: ')
    for k,v in pseudoperson_map.items():
        symbol = HL7.PSEUDOPERSON_SYMBOL_MAP[k]
        count  = len(v)
        print('\t{0} => {1}'.format(symbol, count))
    print()
    
    # compute date tuples; the file data objects contain datetimes,
    # so keep only the date part    
    for pp_category, obj_list in pseudoperson_map.items():
        tuple_dist = defaultdict(int)
        for obj in obj_list:
            # extract all tuple dates from this record
            dt_list = []
            for f in HL7.DATE_TUPLE_FIELDS:
                d = getattr(obj, f, None)
                if d is not None and len(d) > 0:
                    # convert the YYYY-mm-dd part to a datetime
                    dt = datetime.datetime.strptime(d[:10], '%Y-%m-%d')
                    dt_list.append(dt)

            # first_elec_submit_dt is always present, guaranteed by preprocessor
            assert len(dt_list) > 0

            # compute the tuple anchor date
            if _TUPLE_ANCHOR_EARLIEST == _TUPLE_ANCHOR:
                anchor_dt = min(dt_list)
            else:
                anchor_dt = _median_datetime(dt_list)

            # build tuple as day diffs relative to the anchor date
            tup = []
            dt_list_index = 0
            for f in HL7.DATE_TUPLE_FIELDS:
                d = getattr(obj, f, None)
                if d is None or 0 == len(d):
                    tup.append(None)
                else:
                    assert dt_list_index >= 0
                    assert dt_list_index < len(dt_list)
                    dt = dt_list[dt_list_index]
                    dt_list_index += 1
                    day_diff = (dt - anchor_dt).days
                    tup.append(day_diff)
            tup = tuple(tup)

            if _TUPLE_ANCHOR_EARLIEST == _TUPLE_ANCHOR:
                # tuple components should either be None or nonnegative
                for c in tup:
                    if c is not None:
                        assert c >= 0

            # if tup is not None:
            #     print('tup: {0}'.format(tup))
            #     print('\t       report_dt_new: {0}'.format(obj.report_dt_new))
            #     print('\tfirst_elec_submit_dt: {0}'.format(obj.first_elec_submit_dt))
            #     print('\t             diag_dt: {0}'.format(obj.diag_dt))
            #     print('\t             died_dt: {0}'.format(obj.died_dt))
            #     print('\t       hosp_admit_dt: {0}'.format(obj.hosp_admit_dt))
            #     print('\t    illness_onset_dt: {0}'.format(obj.illness_onset_dt))
            #     print('\t     invest_start_dt: {0}'.format(obj.invest_start_dt))

            # update the tuple distribution for this pseudoperson type
            tuple_dist[tup] += 1

        pp = _to_pseudoperson_obj(pp_category, tuple_dist)
        _state.pseudoperson_map[pp_category] = pp

        _print_top_pseudoperson_tuples(pp_category)


###############################################################################
def plot_pseudoperson_distributions(distributions_only = False):
    """
    Plot pseudoperson tuple distributions, ECDFs and inverses.
    """

    for k,pp in _state.pseudoperson_map.items():
        tuple_index_samples = []
        for index, count in pp.tuple_index_dist.items():
            for k in range(count):
                tuple_index_samples.append(index)

        plots_hl7.plot_tuple_pdf_ecdf(len(pp.tuple_index_dist),
                                      tuple_index_samples,
                                      pp.tuple_ecdf,
                                      distributions_only)


###############################################################################
def _merge_files(filepath_list, rng):
    """
    Merge the data in the filepath list and write the merged data to
    a temp file.

    All files in the list are assumed to exist. The front-end code checks each
    file for existence.
    """

    DATE_FORMAT = '%Y-%m-%d'
    header = None
    num_cols = None

    data_dict = defaultdict(list)
    start_date = None
    end_date   = None

    for filepath in filepath_list:
        with open(filepath, 'rt') as infile:
            # load file in single gulp
            raw = infile.read()
            lines = raw.split('\n')
            file_start_date = None
            file_end_date = None
            for i,line in enumerate(lines):
                # skip header line
                if 0 == i:
                    if header is None:
                        # save header
                        header = line.strip()
                        num_cols = len(header.split(','))
                    continue
                text = line.strip()
                if 0 == len(text):
                    continue

                # get components from the next line of data
                item_list = text.split(',')

                str_ref_date = item_list[HL7.OFFSET_REF_DATE]
                assert str_ref_date is not None
                assert len(str_ref_date) > 0
                
                # this could be an empty string
                str_count  = item_list[HL7.OFFSET_COUNT]

                # update the start_date and record the final date seen
                if file_start_date is None:
                    file_start_date = datetime.datetime.strptime(str_ref_date[:10],
                                                                 DATE_FORMAT)
                    if start_date is None or file_start_date < start_date:
                        start_date = file_start_date

                # record the final date seen, which will be used to determine
                # the end date for this file
                file_end_date = str_ref_date

                # if COUNT > 0 keep the line, otherwise discard
                if str_count is not None and len(str_count) > 0:
                    count = float(str_count)                
                    if 0 == count:
                        continue
                    else:
                        data_dict[str_ref_date].append(text)
                else:
                    # assume this line is corrupted
                    continue

            # update the end_date
            end_date_this_file = datetime.datetime.strptime(file_end_date[:10],
                                                            DATE_FORMAT)
            if end_date is None or end_date_this_file > end_date:
                end_date = end_date_this_file

    assert start_date is not None
    assert end_date is not None

    # generate a temporary filename
    hash_string = '|'.join(filepath_list)+str(time.time())
    sha256 = hashlib.sha256()
    sha256.update(str.encode(hash_string))
    outfile_name = sha256.digest().hex()[:32] + '.csv'

    if _TRACE:
        print('From _merge_files(): ')
        print('\t      Start Date: {0}'.format(start_date))
        print('\t        End Date: {0}'.format(end_date))
        print('\tMerged File Name: {0}'.format(outfile_name))

    # write the merged data to disk, using COUNT==0 records for missing data
    cur_date = start_date
    with open(outfile_name, 'wt') as outfile:
        # write the header
        outfile.write('{0}\n'.format(header))
        while cur_date <= end_date:
            str_cur_date = cur_date.strftime(DATE_FORMAT)
            if str_cur_date in data_dict:
                # file data exists for the current date, so write it out
                record_list = data_dict[str_cur_date]
                # scramble the merged records for this date, so that the
                # data blocks from each file do not always appear in the
                # same order
                rng.shuffle(record_list)
                for r in record_list:

                    # verify that the reference date is correct
                    items = r.split(',')
                    str_check_date = items[HL7.OFFSET_REF_DATE][:10]
                    if str_check_date != str_cur_date:
                        print('\n*** File merge problem ***')
                        print('str_check_date: "{0}"'.format(str_check_date))
                        print('  str_cur_date: "{0}"'.format(str_cur_date))
                        for i,item in enumerate(items):
                            print('[{0:2d}]: {1}'.format(i, item))
                        raise SystemExit('_merge_files: invalid check date')

                    outfile.write('{0}\n'.format(r))
            else:
                # write an empty record for the current date
                empty_line_template = HL7.EMPTY_LINE
                outfile.write(empty_line_template.format(str_cur_date))

            # increment to the next calendar day (handles leap years)
            cur_date += datetime.timedelta(days=1)
            
    return outfile_name


###############################################################################
def _empty_record(ref_date):
    """
    Create and return an empty _HL7_Record for the given ref date.
    """
    
    return _HL7_Record(
        event_date            = ref_date,
        count                 = 0,
        age                   = None,
        sex                   = None,
        race                  = None,
        ethnicity             = None,
        case_status           = None,
        county                = None,
        pregnant              = None,
        first_elec_submit_dt  = None,
        diag_dt               = None,
        died_dt               = None,
        report_dt_new         = None,
        hosp_admit_dt         = None,
        illness_onset_dt      = None,
        invest_start_dt       = None,
    )


###############################################################################
def _has_valid_prefix(test_str, hl7_map):
    """
    Test all keys in an HL7 variable map and determine whether 'test_str'
    has a key as a prefix. If so, return the key, otherwise return None.
    
    Searches only for word prefixes.
    """

    # append a space to the test string
    search_str = test_str + ' '

    if _TRACE:
        print('_has_valid_prefix: searching for "{0}" as a prefix: '.
              format(test_str))

    # the keys are the official values
    for k,v in hl7_map.items():
        if k.startswith(search_str):
            # initial words of search_str are a prefix of the official value
            if _TRACE:
               print('\tfound "{0}" as a prefix of "{1}"'.
                     format(test_str, k))
            return k
        elif search_str.startswith(k):
            # the official value is a prefix of search_str
            if _TRACE:
                print('\tfound "{0}" as a prefix of "{1}"'.
                      format(k, search_str))
            return k

    if _TRACE:
        print('\tnot found: "{0}" is not a valid prefix'.format(search_str))
            
    return None


###############################################################################
def _compute_age(str_age, str_age_units, problem_counts):
    """
    Attempt to compute the age in years from the given age and age units.
    Return None for the age value if the calculation cannot be performed.
    """

    age = None

    valid_age = str_age is not None and len(str_age) > 0 and \
        HL7.DATE_UNKNOWN != str_age
    valid_age_units = str_age_units is not None and len(str_age_units) > 0 and \
        str_age_units in HL7.COMPUTABLE_AGE_UNITS
        
    if valid_age and valid_age_units:
        age = float(str_age)
        # convert age to years according to age_units
        if HL7.AGE_UNITS_YEARS == str_age_units:
            # years, no conversion needed
            pass
        elif HL7.AGE_UNITS_MONTHS == str_age_units:
            # months
            age = age / 12.0
        elif HL7.AGE_UNITS_WEEKS == str_age_units:
            # weeks
            age = age / 52.0
        elif HL7.AGE_UNITS_DAYS == str_age_units:
            # days
            age = age / 365.0
        else:
            print('Found unhandled age_units "{0}"'.format(str_age_units))
            problem_counts[_UNHANDLED_AGE_UNITS] += 1
            # override the age value, cannot interpret
            age = None

    return age


###############################################################################
def _compute_age_from_birthdate(str_birthdate, str_ref_date, problem_counts):
    """
    Attempt to compute the age in years from the birthdate and reference date.
    Return None for the age value if the calculation cannot be performed.
    """

    str1_ok = str_birthdate is not None and len(str_birthdate) > 0 and \
        HL7.DATE_UNKNOWN != str_birthdate
    str2_ok = str_ref_date is not None and len(str_ref_date) > 0 and \
        HL7.DATE_UNKNOWN != str_ref_date
    
    age = None    
    if str1_ok and str2_ok:
        # Found a birthdate, so decode it and compute the age in years as the
        # difference between the ref_date and the birthdate. The birthdate
        # format seems to be YYYYmmdd. Have also seeen this format:
        # YYYYmmdd000000.000. Extract the YYYYmmdd portion and decode that.
        str_birthdate = str_birthdate[:8]
        bd = datetime.datetime.strptime(str_birthdate, '%Y%m%d')
        rd = datetime.datetime.strptime(str_ref_date, '%Y-%m-%d')
        if rd < bd:
            # negative age is invalid
            problem_counts[_MISSING_AGE] += 1
            age = None
        else:
            # this is a datetime.timedelta object
            delta = rd - bd
            # express as a difference in days and convert to years
            age = float(delta.days / 365.0)

    return age


###############################################################################
def _to_hl7_record(str_ref_date,
                   count,
                   str_age,
                   str_age_units,
                   str_birthdate,
                   str_sex,
                   str_race,
                   str_ethnicity,
                   str_case_status,
                   str_report_county,
                   str_subj_county,
                   str_pregnant,
                   str_diag_dt,
                   str_died_dt,
                   str_report_dt_new,
                   str_first_elec_submit_dt,
                   str_hosp_admit_dt,
                   str_illness_onset_dt,
                   str_invest_start_dt):
    """
    Decode the data from a single line of an input file.
    """

    age         = None
    sex         = None
    race        = None
    ethnicity   = None
    county      = None
    case_status = None
    pregnant    = None
    problem_counts = defaultdict(int)

    # compute age in years and check if within bounds
    age = _compute_age(str_age, str_age_units, problem_counts)
    if age is None:
        # try to determine from birthdate and ref date
        age = _compute_age_from_birthdate(str_birthdate,
                                          str_ref_date,
                                          problem_counts)        
    if age is None:
        # could not determine age
        age = HL7.AGE_UNKNOWN
    else:
        # an age value was computed, check if within bounds
        assert float == type(age)
        if age < 0.0 or age > HL7.AGE_MAX_FLOAT:
            age = HL7.AGE_UNKNOWN

    # remap sex strings to floats
    if str_sex is None or 0 == len(str_sex):
        problem_counts[_MISSING_SEX] += 1
        sex = HL7.SEX_MAP[HL7.SEX_UNKNOWN]
    else:
        if str_sex not in HL7.SEX_MAP:
            # invalid value
            problem_counts[_INVALID_SEX] += 1
            sex = HL7.SEX_MAP[HL7.SEX_UNKNOWN]
            if _TRACE:
                print('Invalid sex value: "{0}"'.format(str_sex))
        else:
            sex = HL7.SEX_MAP[str_sex]

    # the race is a string that could be a semicolon-concatenated list of
    # individual race identifiers
    # return race as a string
    if str_race is None or 0 == len(str_race):
        problem_counts[_MISSING_RACE] += 1
        race = str(HL7.RACE_UNKNOWN)
    else:
        # needs remapping later
        race = str_race

    # remap ethnicity strings to floats
    if str_ethnicity is None or 0 == len(str_ethnicity):
        problem_counts[_MISSING_ETHNICITY] += 1
        ethnicity = HL7.ETHNICITY_MAP[HL7.ETHNICITY_UNKNOWN]
    else:
        if str_ethnicity not in HL7.ETHNICITY_MAP:
            # check to see if str_ethnicity starts with a known value
            key = _has_valid_prefix(str_ethnicity, HL7.ETHNICITY_MAP)
            if key is not None:
                ethnicity = HL7.ETHNICITY_MAP[key]
            else:
                # invalid value
                problem_counts[_INVALID_ETHNICITY] += 1
                ethnicity = HL7.ETHNICITY_MAP[HL7.ETHNICITY_UNKNOWN]
                if _TRACE:
                    print('Invalid ethnicity value: "{0}"'.format(str_ethnicity))
        else:
            ethnicity = HL7.ETHNICITY_MAP[str_ethnicity]

    # remap case status strings to floats
    if str_case_status is None or 0 == len(str_case_status):
        problem_counts[_MISSING_CASE_STATUS] += 1
        case_status = HL7.CASE_STATUS_MAP[HL7.CASE_STATUS_UNKNOWN]
    else:
        if str_case_status not in HL7.CASE_STATUS_MAP:
            # check to see if str_case_status starts with a known value
            key = _has_valid_prefix(str_case_status, HL7.CASE_STATUS_MAP)
            if key is not None:
                case_status = HL7.CASE_STATUS_MAP[key]
            else:
                # invalid value
                problem_counts[_INVALID_CASE_STATUS] += 1
                case_status = HL7.CASE_STATUS_MAP[HL7.CASE_STATUS_UNKNOWN]
                if _TRACE:
                    print('Invalid case status value: "{0}"'.format(str_case_status))
        else:
            case_status = HL7.CASE_STATUS_MAP[str_case_status]

    # use str_subj_county if present, otherwise use str_report_county
    # return the county as a STRING, not a float
    bad_subj_county = str_subj_county is None or 0 == len(str_subj_county)
    bad_report_county = str_report_county is None or 0 == len(str_report_county)
    if bad_subj_county and bad_report_county:
        problem_counts[_MISSING_COUNTY] += 1
        county = str(HL7.COUNTY_UNKNOWN)
    elif not bad_subj_county:
        county = str_subj_county
    else:
        # needs remapping later
        county = str_report_county

    # remap pregnant strings to floats
    if str_pregnant is None or 0 == len(str_pregnant):
        problem_counts[_MISSING_PREGNANT] += 1
        pregnant = HL7.PREGNANT_MAP[HL7.PREGNANT_UNKNOWN]
    else:
        if str_pregnant not in HL7.PREGNANT_MAP:
            # check to see if str_pregnant starts with a known value
            key = _has_valid_prefix(str_pregnant, HL7.PREGNANT_MAP)
            if key is not None:
                pregnant = HL7.PREGNANT_MAP[key]
            else:
                # invalid value
                problem_counts[_INVALID_PREGNANT] += 1
                pregnant = HL7.PREGNANT_MAP[HL7.PREGNANT_UNKNOWN]
                if _TRACE:
                    print('Invalid pregnant value: "{0}"'.format(str_pregnant))
        else:
            pregnant = HL7.PREGNANT_MAP[str_pregnant]
        
    assert str_ref_date is not None
    assert len(str_ref_date) > 0
    assert count is not None
    assert age is not None
    assert sex is not None
    assert race is not None
    assert ethnicity is not None
    assert case_status is not None
    assert county is not None
    assert pregnant is not None

    # check categorical variables with a finite set of options
    assert sex in HL7.SEX_MAP_VALUES
    assert ethnicity in HL7.ETHNICITY_MAP_VALUES
    assert case_status in HL7.CASE_STATUS_MAP_VALUES
    assert pregnant in HL7.PREGNANT_MAP_VALUES

    assert isinstance(race, str)
    assert isinstance(county, str)

    record = _HL7_Record(
        event_date           = str_ref_date,
        count                = count,
        age                  = float(age),
        sex                  = float(sex),
        race                 = race,
        ethnicity            = float(ethnicity),
        case_status          = float(case_status),
        county               = county,
        pregnant             = float(pregnant),
        diag_dt              = str_diag_dt,
        died_dt              = str_died_dt,
        first_elec_submit_dt = str_first_elec_submit_dt,
        report_dt_new        = str_report_dt_new,
        hosp_admit_dt        = str_hosp_admit_dt,
        illness_onset_dt     = str_illness_onset_dt,
        invest_start_dt      = str_invest_start_dt
    )

    return record, problem_counts


###############################################################################
def _occurs_after(dt1, str_date2):
    """
    Convert str_date2 to a datetime and return True if it occurs after date1,
    which is already a datetime.
    """

    if str_date2 is not None and len(str_date2) > 0:
        dt2 = datetime.datetime.strptime(str_date2[:10], '%Y-%m-%d')
        if dt2 > dt1:
            return True

    return False


###############################################################################
def _occurs_before(dt1, str_date2):
    """
    Convert str_date2 to a datetime and return True if it occurs before date1,
    which is already a datetime.
    """

    if str_date2 is not None and len(str_date2) > 0:
        dt2 = datetime.datetime.strptime(str_date2[:10], '%Y-%m-%d')
        if dt2 < dt1:
            return True

    return False


###############################################################################
def _load_single_file(filepath):
    """
    Load a single preprocessed HL7 data file.
    """

    file_data = []
    line_count = 0
    lines_with_data = 0
    zero_count_fields = 0
    bad_record_count = 0
    count_sum = 0
    county_set = set()

    total_problem_counts = defaultdict(int)
    total_problem_counts[_MISSING_AGE] = 0
    total_problem_counts[_MISSING_AGE_UNITS] = 0
    total_problem_counts[_MISSING_SEX] = 0
    total_problem_counts[_MISSING_RACE] = 0
    total_problem_counts[_MISSING_ETHNICITY] = 0
    total_problem_counts[_MISSING_COUNTY] = 0
    total_problem_counts[_MISSING_CASE_STATUS] = 0
    total_problem_counts[_MISSING_PREGNANT] = 0
    total_problem_counts[_INVALID_SEX] = 0
    total_problem_counts[_INVALID_RACE] = 0
    total_problem_counts[_INVALID_ETHNICITY] = 0
    total_problem_counts[_INVALID_CASE_STATUS] = 0
    total_problem_counts[_INVALID_PREGNANT] = 0
    total_problem_counts[_UNHANDLED_AGE_UNITS] = 0
    total_problem_counts[_NOT_A_CASE] = 0
    total_problem_counts[_RESULTS_NOT_OBTAINED] = 0

    records = []

    if _TRACE:
        race_values = set()
        county_values = set()
        pregnant_values = set()
        nrs_dict = defaultdict(int)

    is_json_file = False
    if filepath.lower().endswith('csv'):

        with open(filepath, 'rt') as infile:
            # load file in single gulp
            raw = infile.read()
            lines = raw.split('\n')

            # save and skip the header line
            _state.preprocessed_header_strings = lines[0].lower().split(',')
            lines = lines[1:]

            # strip newlines
            lines = [line.strip() for line in lines]

            # prune empty lines
            records = [line for line in lines if len(line) > 0]

    elif filepath.lower().endswith('json'):
        
        is_json_file = True
        with open(filepath, 'rt') as infile:
            raw = infile.read()
            records = json.loads(raw)

    bad_record_count  = 0
    zero_count_fields = 0

    for obj in records:
        if is_json_file:
            str_report_dt_new        = None
            str_count                = None
            str_age                  = None
            str_age_units            = None
            str_sex                  = None
            str_ethnicity            = None
            str_race_mapped          = None
            str_case_status          = None
            str_birthdate            = None
            str_notif_result_status  = None
            str_pregnant             = None
            str_report_county        = None
            str_first_elec_submit_dt = None
            str_subj_county          = None
            str_diag_dt              = None
            str_died_dt              = None
            str_earliest_cnty_dt     = None
            str_earliest_state_dt    = None
            str_hosp_admit_dt        = None
            str_illness_onset_dt     = None
            str_invest_start_dt      = None
            str_phd_notif_dt         = None

            if HL7.FIELD_REPORT_DATE in obj:
                str_report_dt_new = str(obj[HL7.FIELD_REPORT_DATE])
            if HL7.FIELD_COUNT in obj:
                str_count = str(obj[HL7.FIELD_COUNT])
            if HL7.FIELD_AGE in obj:
                str_age = str(obj[HL7.FIELD_AGE])
            if HL7.FIELD_AGE_UNITS in obj:
                str_age_units = str(obj[HL7.FIELD_AGE_UNITS])
            if HL7.FIELD_SEX in obj:
                str_sex = str(obj[HL7.FIELD_SEX])
            if HL7.FIELD_ETHNICITY in obj:
                str_ethnicity = str(obj[HL7.FIELD_ETHNICITY])
            if HL7.FIELD_RACE in obj:
                str_race_mapped = str(obj[HL7.FIELD_RACE])
            if HL7.FIELD_CASE_STATUS in obj:
                str_case_status = str(obj[HL7.FIELD_CASE_STATUS])
            if HL7.FIELD_BIRTH_DATE in obj:
                str_birthdate = str(obj[HL7.FIELD_BIRTH_DATE])
            if HL7.FIELD_NRS in obj:
                str_notif_result_status = str(obj[HL7.FIELD_NRS])
            if HL7.FIELD_PREGNANT in obj:
                str_pregnant = str(obj[HL7.FIELD_PREGNANT])
            if HL7.FIELD_COUNTY in obj:
                str_report_county = str(obj[HL7.FIELD_COUNTY])
            if HL7.FIELD_DATE1 in obj:
                str_first_elec_submit_dt = str(obj[HL7.FIELD_DATE1])
            if HL7.FIELD_SUBJ_COUNTY in obj:
                str_subj_county = str(obj[HL7.FIELD_SUBJ_COUNTY])
            if HL7.FIELD_DATE2 in obj:
                str_diag_dt = str(obj[HL7.FIELD_DATE2])
            if HL7.FIELD_DATE3 in obj:
                str_died_dt = str(obj[HL7.FIELD_DATE3])
            if HL7.FIELD_DATE4 in obj:
                str_earliest_cnty_dt = str(obj[HL7.FIELD_DATE4])
            if HL7.FIELD_DATE5 in obj:
                str_earliest_state_dt = str(obj[HL7.FIELD_DATE5])
            if HL7.FIELD_DATE6 in obj:
                str_hosp_admit_dt = str(obj[HL7.FIELD_DATE6])
            if HL7.FIELD_DATE7 in obj:
                str_illness_onset_dt = str(obj[HL7.FIELD_DATE7])
            if HL7.FIELD_DATE8 in obj:
                str_invest_start_dt = str(obj[HL7.FIELD_DATE8])
            if HL7.FIELD_DATE9 in obj:
                str_phd_notif_dt = str(obj[HL7.FIELD_DATE9])

        else:
            # obj is a string of comma-separated values
            vars = obj.split(',')

            # strip extraneous whitespace and convert to lowercase
            for i in range(len(vars)):
                if vars[i] is not None and len(vars[i]) > 0:
                    vars[i] = vars[i].strip().lower()
            
            # [0]  Report_Dt_new
            # [1]  COUNT
            # [2]  Age
            # [3]  Age_Units
            # [4]  Sex
            # [5]  Ethnicity_txt
            # [6]  race_mapped
            # [7]  Case_Status_txt
            # [8]  Birth_Date_str
            # [9]  Notif_Result_Status
            # [10] Pregnant
            # [11] Report_County
            # [12] First_Elec_Submit_Dt
            # [13] Subj_County
            # [14] Diag_Dt
            # [15] Died_Dt
            # [16] Earliest_Cnty_Dt
            # [17] Earliest_State_Dt
            # [18] Hosp_Admit_Dt
            # [19] Illness_Onset_Dt
            # [20] Invest_Start_Dt
            # [21] PHD_Notif_Dt

            str_report_dt_new        = vars[0]
            str_count                = vars[1]
            str_age                  = vars[2]
            str_age_units            = vars[3]
            str_sex                  = vars[4]
            str_ethnicity            = vars[5]
            str_race_mapped          = vars[6]
            str_case_status          = vars[7]
            str_birthdate            = vars[8]
            str_notif_result_status  = vars[9]
            str_pregnant             = vars[10]
            str_report_county        = vars[11]
            str_first_elec_submit_dt = vars[12]
            str_subj_county          = vars[13]
            str_diag_dt              = vars[14]
            str_died_dt              = vars[15]
            str_earliest_cnty_dt     = vars[16]
            str_earliest_state_dt    = vars[17]
            str_hosp_admit_dt        = vars[18]
            str_illness_onset_dt     = vars[19]
            str_invest_start_dt      = vars[20]
            str_phd_notif_dt         = vars[21]

        # Keep the YYYY-MM-DD portion of the reference date for the
        # signal processing code. The preprocessor ensures that this
        # date is always present.
        if HL7.FIELD_DATE1 == HL7.FIELD_REF_DATE:
            str_ref_date = str_first_elec_submit_dt[:10]
        else:
            str_ref_date = str_report_dt_new[:10]

        if str_ref_date is None or 0 == len(str_ref_date):
            if HL7.FIELD_DATE1 == HL7.FIELD_REF_DATE:
                print('REF DATE IS FIELD_DATE1')
            else:
                print('REF DATE IS report_dt_new')
            print(obj)
            
        assert str_ref_date is not None
        assert len(str_ref_date) > 0

        # zero count records have no data
        if str_count is not None and len(str_count) > 0:
            count = float(str_count)
            if 0 == count:
                zero_count_fields += 1
                file_data.append(_empty_record(str_ref_date))
                continue
        else:
            # assume this line is corrupted
            bad_record_count += 1
            file_data.append(_empty_record(str_ref_date))
            continue

        if _TRACE:
            race_values.add(str_race_mapped)
            county_values.add('{0};{1}'.format(str_subj_county, str_report_county))
            pregnant_values.add(str_pregnant)
            nrs_dict[str_notif_result_status] += 1
        
        # ignore any records having either of these conditions:
        #     1. str_case_status == 'not a case'
        #     2. str_notif_result_status == 'x'
        if HL7.CASE_STATUS_NOT == str_case_status:
            total_problem_counts[_NOT_A_CASE] += 1
            file_data.append(_empty_record(str_ref_date))
            continue
        elif HL7.NOTIF_RESULT_STATUS_NOT_OBTAINED == str_notif_result_status:
            total_problem_counts[_RESULTS_NOT_OBTAINED] += 1
            file_data.append(_empty_record(str_ref_date))
            continue

        # ignore any dates prior to _MIN_DATE
        if _occurs_before(_MIN_DATE, str_report_dt_new):
            str_report_dt_new = _EMPTY_STRING
        if _occurs_before(_MIN_DATE, str_diag_dt):
            str_diag_dt = _EMPTY_STRING
        if _occurs_before(_MIN_DATE, str_died_dt):
            str_died_dt = _EMPTY_STRING
        if _occurs_before(_MIN_DATE, str_report_dt_new):
            str_report_dt_new = _EMPTY_STRING
        if _occurs_before(_MIN_DATE, str_hosp_admit_dt):
            str_hosp_admit_dt = _EMPTY_STRING
        if _occurs_before(_MIN_DATE, str_illness_onset_dt):
            str_illness_onset_dt = _EMPTY_STRING
        if _occurs_before(_MIN_DATE, str_invest_start_dt):
            str_invest_start_dt = _EMPTY_STRING

        # if _occurs_after(_MAX_DATE, str_report_dt_new):
        #     str_report_dt_new = _EMPTY_STRING
        # if _occurs_after(_MAX_DATE, str_diag_dt):
        #     str_diag_dt = _EMPTY_STRING
        # if _occurs_after(_MAX_DATE, str_died_dt):
        #     str_died_dt = _EMPTY_STRING
        # if _occurs_after(_MAX_DATE, str_report_dt_new):
        #     str_report_dt_new = _EMPTY_STRING
        # if _occurs_after(_MAX_DATE, str_hosp_admit_dt):
        #     str_hosp_admit_dt = _EMPTY_STRING
        # if _occurs_after(_MAX_DATE, str_illness_onset_dt):
        #     str_illness_onset_dt = _EMPTY_STRING
        # if _occurs_after(_MAX_DATE, str_invest_start_dt):
        #     str_invest_start_dt = _EMPTY_STRING

        # IMPORTANT: load all data BEFORE jittering any dates
        #            the age is calculated as birth_date - ref_date
        #            if no age value is provided; altering the ref_date
        #            would therefore alter the age distribution

        hl7_record, problem_counts = _to_hl7_record(str_ref_date,
                                                    count,
                                                    str_age,
                                                    str_age_units,
                                                    str_birthdate,
                                                    str_sex,
                                                    str_race_mapped,
                                                    str_ethnicity,
                                                    str_case_status,
                                                    str_report_county,
                                                    str_subj_county,
                                                    str_pregnant,
                                                    str_diag_dt,
                                                    str_died_dt,
                                                    str_report_dt_new,
                                                    str_first_elec_submit_dt,
                                                    str_hosp_admit_dt,
                                                    str_illness_onset_dt,
                                                    str_invest_start_dt)


        for k,v in problem_counts.items():
            total_problem_counts[k] += v

        # count the number of different county values
        if hl7_record.county is not None:
            county_set.add(hl7_record.county)

        # update some stats
        lines_with_data += 1
        count_sum += int(hl7_record.count)

        file_data.append(hl7_record)

    # required width for line counts
    line_count = len(records)
    width = len(str(line_count))    

    if _TRACE:
        print('Race values: ')
        for v in race_values:
            print('\t"{0}"'.format(v))
        print()

        print('County values (subj_county;report_county): ')
        for v in county_values:
            print('\t"{0}"'.format(v))
        print()

        print('Pregnant values: ')
        for v in pregnant_values:
            print('\t"{0}"'.format(v))
        print()

        print('notif_result_status values and counts: ')
        for k,v in nrs_dict.items():
            print('\t{0} => {1}'.format(k, v))
        print()

    print('Information for file "{0}":'.format(filepath))
    print('\t      Start Elect Submit Date: {0}'.format(file_data[0].event_date))
    print('\t        End Elect Submit Date: {0}'.format(file_data[-1].event_date))
    print('\tLine count (excluding header): {0}'.format(line_count, width))
    print('\t  Lines with complete records: {0:{1}d}'.format(lines_with_data, width))
    print('\t      Sum of all count fields: {0:{1}d}'.format(count_sum, width))
    print('\t        Lines with COUNT == 0: {0:{1}d}'.format(zero_count_fields, width))
    print('\t           NOT_A_CASE entries: {0:{1}d}'.format(total_problem_counts[_NOT_A_CASE],
                                                             width))
    print('\t RESULTS_NOT_OBTAINED entries: {0:{1}d}'.format(total_problem_counts[_RESULTS_NOT_OBTAINED],
                                                             width))
    print('\t    No. of different counties: {0:{1}d}'.format(len(county_set), width))
    print('\t                  Bad records: {0:{1}d}'.format(bad_record_count, width))
    print('\t          Missing AGE entries: {0:{1}d}'.format(total_problem_counts[_MISSING_AGE],
                                                             width))
    print('\t    Missing AGE_UNITS entries: {0:{1}d}'.format(total_problem_counts[_MISSING_AGE_UNITS],
                                                             width))
    print('\t  Unhandled AGE_UNITS entries: {0:{1}d}'.format(total_problem_counts[_UNHANDLED_AGE_UNITS],
                                                             width))
    print('\t          Missing SEX entries: {0:{1}d}'.format(total_problem_counts[_MISSING_SEX],
                                                             width))
    print('\t         Missing RACE entries: {0:{1}d}'.format(total_problem_counts[_MISSING_RACE],
                                                             width))
    print('\t    Missing ETHNICITY entries: {0:{1}d}'.format(total_problem_counts[_MISSING_ETHNICITY],
                                                             width))
    print('\t       Missing COUNTY entries: {0:{1}d}'.format(total_problem_counts[_MISSING_COUNTY],
                                                             width))
    print('\t  Missing CASE_STATUS entries: {0:{1}d}'.format(total_problem_counts[_MISSING_CASE_STATUS],
                                                             width))
    print('\t     Missing PREGNANT entries: {0:{1}d}'.format(total_problem_counts[_MISSING_PREGNANT],
                                                             width))
    print('\t          Invalid SEX entries: {0:{1}d}'.format(total_problem_counts[_INVALID_SEX],
                                                             width))
    print('\t         Invalid RACE entries: {0:{1}d}'.format(total_problem_counts[_INVALID_RACE],
                                                             width))
    print('\t    Invalid ETHNICITY entries: {0:{1}d}'.format(total_problem_counts[_INVALID_ETHNICITY],
                                                             width))
    print('\t  Invalid CASE_STATUS entries: {0:{1}d}'.format(total_problem_counts[_INVALID_CASE_STATUS],
                                                             width))
    print('\t     Invalid PREGNANT entries: {0:{1}d}'.format(total_problem_counts[_INVALID_PREGNANT],
                                                             width))

    # the generator cannot handle datasets with only a single valid record
    if lines_with_data <= 1:
        print('\n*** Insufficient source data ***')
        if 0 == lines_with_data:
            record_word = 'records'
        else:
            record_word = 'record'
        print('\tSource dataset contains {0} complete {1}.'.
              format(lines_with_data, record_word))
        return None, None

    return file_data, lines_with_data


###############################################################################
def _to_pdf(bins):
    """
    Normalize a list of bin counts to get an estimated prob distribution.
    """
    
    the_sum = 0.0
    for count in bins:
        the_sum += count
    pdf = np.zeros(len(bins))
    for i in range(len(bins)):
        pdf[i] = bins[i] / the_sum
    return pdf


###############################################################################
def _zero_fill_missing_bins(min_val, max_val, hist_dict):
    """
    For the histogram (stored as a dict) in 'hist_dict', add zero-value bins
    at each integer in the range [min_val, max_val] if those bins are not
    already present. This is to get a contiguous range of bins for the
    histogram.
    """
    
    for i in range(min_val, max_val+1):
        if i not in hist_dict:
            hist_dict[float(i)] = 0


###############################################################################
def _do_consistency_checks():
    """
    Check the data arrays to ensure that all arrays have the expected number
    of elements. Also check the PDFs.
    """

    age_counts = 0
    sex_counts = 0
    race_counts = 0
    ethnicity_counts = 0
    county_counts = 0
    case_status_counts = 0
    pregnant_counts = 0

    for k in sorted(_state.age_dict.keys()):
        age_counts += _state.age_dict[k]
    for k in sorted(_state.sex_dict.keys()):
        sex_counts += _state.sex_dict[k]
    for k in sorted(_state.race_dict.keys()):
        race_counts += _state.race_dict[k]
    for k in sorted(_state.ethnicity_dict.keys()):
        ethnicity_counts += _state.ethnicity_dict[k]
    for k in sorted(_state.case_status_dict.keys()):
        case_status_counts += _state.case_status_dict[k]
    for k in sorted(_state.county_dict.keys()):
        county_counts += _state.county_dict[k]    
    for k in sorted(_state.pregnant_dict.keys()):
        pregnant_counts += _state.pregnant_dict[k]

    if _TRACE:
        print()
        print('AGE: ')
        for i, k in enumerate(sorted(_state.age_dict.keys())):
            print('\t[{0:3d}]: {1} => {2}'.format(i, k, _state.age_dict[k]))
        print('SEX: ')
        for i, k in enumerate(sorted(_state.sex_dict.keys())):
            print('\t[{0}]: {1} => {2}'.format(i, k, _state.sex_dict[k]))
        print('RACE: ')
        for i, k in enumerate(sorted(_state.race_dict.keys())):
            print('\t[{0:3d}]: {1} => {2}'.format(i, k, _state.race_dict[k]))
        print('ETHNICITY: ')
        for i, k in enumerate(sorted(_state.ethnicity_dict.keys())):
            print('\t[{0}]: {1} => {2}'.format(i, k, _state.ethnicity_dict[k]))
        print('CASE STATUS: ')
        for i, k in enumerate(sorted(_state.case_status_dict.keys())):
            print('\t[{0}]: {1} => {2}'.format(i, k, _state.case_status_dict[k]))
        print('COUNTY: ')
        for i, k in enumerate(sorted(_state.county_dict.keys())):
            print('\t[{0:3d}]: {1} => {2}'.format(i, k, _state.county_dict[k]))
        print('PREGNANT: ')
        for i, k in enumerate(sorted(_state.pregnant_dict.keys())):
            print('\t[{0}]: {1} => {2}'.format(i, k, _state.pregnant_dict[k]))

        print()
        print('                     race_max: {0}'.format(_state.race_max))
        print('                   county_max: {0}'.format(_state.county_max))

        print()
        print('        Nonempty record count: {0}'.format(_state.nonempty_record_count))
        print('        Sum of all sex counts: {0}'.format(sex_counts))
        print('       Sum of all race counts: {0}'.format(race_counts))
        print('  Sum of all ethnicity counts: {0}'.format(ethnicity_counts))
        print('     Sum of all county counts: {0}'.format(county_counts))
        print('Sum of all case status counts: {0}'.format(case_status_counts))
        print('   Sum of all pregnant counts: {0}'.format(pregnant_counts))
        print()

    # compute marginal pdfs from this data
    if age_counts > 0:
        age_pdf = _to_pdf(list(_state.age_dict.values()))
        assert _state.nonempty_record_count == age_counts
    if sex_counts > 0:
        sex_pdf = _to_pdf(list(_state.sex_dict.values()))
        assert _state.nonempty_record_count == sex_counts
    if race_counts > 0:
        race_pdf = _to_pdf(list(_state.race_dict.values()))
        assert _state.nonempty_record_count == race_counts
    if ethnicity_counts > 0:
        ethnicity_pdf = _to_pdf(list(_state.ethnicity_dict.values()))
        assert _state.nonempty_record_count == ethnicity_counts
    if case_status_counts > 0:
        case_status_pdf = _to_pdf(list(_state.case_status_dict.values()))
        assert _state.nonempty_record_count == case_status_counts
    if county_counts > 0:
        county_pdf = _to_pdf(list(_state.county_dict.values()))
        assert _state.nonempty_record_count == county_counts
    if pregnant_counts > 0:
        pregnant_pdf = _to_pdf(list(_state.pregnant_dict.values()))
        assert _state.nonempty_record_count == pregnant_counts

    if _TRACE:
        # check to see that the PDF bin values all sum to 1
        # (bins are contiguous with width == 1)
        sum_age_pdf = np.sum(age_pdf)
        print('              age_pdf bin sum: {0:.3f}'.format(sum_age_pdf))
        sum_sex_pdf = np.sum(sex_pdf)
        print('              sex_pdf bin sum: {0:.3f}'.format(sum_sex_pdf))
        sum_race_pdf = np.sum(race_pdf)
        print('             race_pdf bin sum: {0:.3f}'.format(sum_race_pdf))
        sum_ethnicity_pdf = np.sum(ethnicity_pdf)
        print('        ethnicity_pdf_bin_sum: {0:.3f}'.format(sum_ethnicity_pdf))
        sum_case_status_pdf = np.sum(case_status_pdf)
        print('      case_status_pdf_bin_sum: {0:.3f}'.format(sum_case_status_pdf))
        sum_county_pdf = np.sum(county_pdf)
        print('           county_pdf_bin_sum: {0:.3f}'.format(sum_county_pdf))
        sum_pregnant_pdf = np.sum(pregnant_pdf)
        print('         pregnant_pdf_bin_sum: {0:.3f}'.format(sum_pregnant_pdf))
        print()


###############################################################################
def _to_sequential_ints(var_data, var_dict):
    """
    """

    var_map = {}
    var_inv_map = {}

    if len(var_data) > 0:
        new_dict = {}

        # take the keys in sorted order, 0 is now the unknown value
        keys = sorted(var_dict.keys())

        index = 0
        for k in keys:
            # assign next key to the next index
            var_map[k] = index
            var_inv_map[index] = k
            new_dict[index] = var_dict[k]
            index += 1

        var_max = len(var_dict)

        # remap the data to use the new values
        new_data = []
        for orig_k in var_data:
            new_k = var_map[orig_k]
            new_data.append(new_k)
        assert len(new_data) == len(var_data)
        
    return new_data, new_dict, var_map, var_inv_map, var_max


###############################################################################
def _remap_data(file_data):
    """
    """

    # number of data entries with the count field value > 0
    _state.nonempty_record_count = 0

    for r in file_data:

        if r.count is None or 0 == r.count:
            continue

        # the file loader ensures that any line with a nonzero count has data
        _state.nonempty_record_count += 1
    
        # process age data
        if r.age is not None:
            if HL7.AGE_UNKNOWN == r.age:
                the_age = HL7.AGE_UNKNOWN_REMAP
            elif r.age > 0.0 and r.age < 1.0:
                # put fractions in the 0 bin
                the_age = 0.0
            elif r.age > HL7.AGE_MAX_FLOAT:
                # any values greater than this go in the unknown bin
                the_age = HL7.AGE_UNKNOWN_REMAP
            else:
                the_age = r.age

            # truncate to the left edge of the age bin
            the_age = float(int(the_age))

            _state.age_data.append(the_age)
            _state.age_dict[the_age] += 1

        # process sex data (no remap needed)
        if r.sex is not None:
            assert r.sex >= HL7.SEX_MIN_FLOAT and r.sex <= HL7.SEX_MAX_FLOAT
            _state.sex_data.append(r.sex)
            _state.sex_dict[r.sex] += 1

        # process race data
        if r.race is not None:
            if HL7.RACE_UNKNOWN == r.race:
                the_race = HL7.RACE_UNKNOWN_REMAP_STR
            else:
                the_race = r.race
            _state.race_data.append(the_race)
            _state.race_dict[the_race] += 1

        # process ethnicity data (no remap needed)
        if r.ethnicity is not None:
            assert r.ethnicity >= HL7.ETHNICITY_MIN_FLOAT and \
                r.ethnicity <= HL7.ETHNICITY_MAX_FLOAT
            _state.ethnicity_data.append(r.ethnicity)
            _state.ethnicity_dict[r.ethnicity] += 1

        # process case_status data (no remap needed)
        if r.case_status is not None:
            assert r.case_status >= HL7.CASE_STATUS_MIN_FLOAT and \
                r.case_status <= HL7.CASE_STATUS_MAX_FLOAT
            _state.case_status_data.append(r.case_status)
            _state.case_status_dict[r.case_status] += 1

        # process county data
        if r.county is not None:
            if HL7.COUNTY_UNKNOWN == r.county:
                the_county = HL7.COUNTY_UNKNOWN_REMAP_STR
            else:
                the_county = r.county
            _state.county_data.append(the_county)
            _state.county_dict[the_county] += 1

        # process pregnant data
        if r.pregnant is not None:
            assert r.pregnant >= HL7.PREGNANT_MIN and \
                r.pregnant <= HL7.PREGNANT_MAX
            _state.pregnant_data.append(r.pregnant)
            _state.pregnant_dict[r.pregnant] += 1

    _zero_fill_missing_bins(HL7.AGE_MIN, HL7.AGE_MAX, _state.age_dict)
    _zero_fill_missing_bins(HL7.SEX_MIN, HL7.SEX_MAX, _state.sex_dict)
    # race handled separately below
    _zero_fill_missing_bins(HL7.ETHNICITY_MIN, HL7.ETHNICITY_MAX, _state.ethnicity_dict)
    _zero_fill_missing_bins(HL7.CASE_STATUS_MIN, HL7.CASE_STATUS_MAX, _state.case_status_dict)
    # county handled separately below
    _zero_fill_missing_bins(HL7.PREGNANT_MIN, HL7.PREGNANT_MAX, _state.pregnant_dict)

    # remap the race and county strings to sequential integers

    _state.race_data, _state.race_dict, _state.race_map, \
        _state.inv_race_map, _state.race_max = _to_sequential_ints(_state.race_data,
                                                                   _state.race_dict)

    _state.county_data, _state.county_dict, _state.county_map, \
        _state.inv_county_map, _state.county_max = _to_sequential_ints(_state.county_data,
                                                                       _state.county_dict)

    _do_consistency_checks()

    if _TRACE:
        # print remapped data
        print()
        print('Remapped RACE data: ')
        inv_race = [(v,k) for k,v in _state.race_map.items()]
        sorted_values = sorted(inv_race, key=lambda x: x[0])
        for v,k in sorted_values:
            print('\t[{0:2}]:\t{1}'.format(v, k))
        print()

        print('Remapped COUNTY data: ')
        inv_county = [(v,k) for k,v in _state.county_map.items()]
        sorted_values = sorted(inv_county, key=lambda x: x[0])
        for v,k in sorted_values:
            print('\t[{0:3}]:\t{1}'.format(v, k))
        print()


###############################################################################
def _compute_cdfs(variable_names, record_count):
    """
    Compute the cumulative distribution functions for the data.
    """

    # use kernel density estimation for the age variable for smaller datasets
    if record_count < 2000:
        _state.ecdf_age = KdeECDF(_state.age_data,
                                  l=HL7.AGE_MIN,
                                  r=HL7.AGE_MAX,
                                  kernel='hybrid')
    else:
        _state.ecdf_age = EmpiricalCDF(_state.age_data)

    _state.ecdf_sex = EmpiricalCDF(_state.sex_data)
    _state.ecdf_race = EmpiricalCDF(_state.race_data)
    _state.ecdf_ethnicity = EmpiricalCDF(_state.ethnicity_data)
    _state.ecdf_case_status = EmpiricalCDF(_state.case_status_data)
    _state.ecdf_county = EmpiricalCDF(_state.county_data)
    _state.ecdf_pregnant = EmpiricalCDF(_state.pregnant_data)

    ecdf_list = []

    # the order of the ECDF objects MUST match the order of the variable_names
    for name in variable_names:
        if HL7.NAME_AGE == name:
            ecdf_list.append(_state.ecdf_age)
        elif HL7.NAME_SEX == name:
            ecdf_list.append(_state.ecdf_sex)
        elif HL7.NAME_RACE == name:
            ecdf_list.append(_state.ecdf_race)
        elif HL7.NAME_ETHNICITY == name:
            ecdf_list.append(_state.ecdf_ethnicity)
        elif HL7.NAME_CASE_STATUS == name:
            ecdf_list.append(_state.ecdf_case_status)
        elif HL7.NAME_COUNTY == name:
            ecdf_list.append(_state.ecdf_county)
        elif HL7.NAME_PREGNANT == name:
            ecdf_list.append(_state.ecdf_pregnant)

    return ecdf_list


###############################################################################
def _validate_names(variable_names):

    for name in variable_names:
        if name not in HL7.NAME_STRING_SET:
            print('*** ERROR ***: invalid variable name "{0)"'.format(name))
            return False

    return True


###############################################################################
def _build_data_list(variable_names):
    """
    """

    data_list = []
    
    for name in variable_names:
        if HL7.NAME_AGE == name:
            data_list.append(_state.age_data)
        elif HL7.NAME_SEX == name:
            data_list.append(_state.sex_data)
        elif HL7.NAME_RACE == name:
            data_list.append(_state.race_data)
        elif HL7.NAME_ETHNICITY == name:
            data_list.append(_state.ethnicity_data)
        elif HL7.NAME_CASE_STATUS == name:
            data_list.append(_state.case_status_data)
        elif HL7.NAME_COUNTY == name:
            data_list.append(_state.county_data)
        elif HL7.NAME_PREGNANT == name:
            data_list.append(_state.pregnant_data)

    return data_list


###############################################################################
def init_model_data(filepath_list, variable_names_in, rng=None):
    """
    Merge and load all input files, build all data structures, and initialize
    everything for a new run.
    """

    ERROR_TUPLE = (None, None, None, None)

    # clear all variables to allow multiple runs
    _state.reset()

    # check the names of the model variables
    if not _validate_names(variable_names_in):
        return ERROR_TUPLE
    variable_names = variable_names_in

    # load or load and merge data files
    if 1 == len(filepath_list):
        file_data, record_count = _load_single_file(filepath_list[0])
    else:
        # Merge the data files and write the merged data to a temp file.
        # Load the temp file to build the data arrays, then delete the
        # temp file.
        assert rng is not None
        temp_filename = _merge_files(filepath_list, rng)
        file_data, record_count = _load_single_file(temp_filename)
        if not _TRACE:
            try:
                os.remove(temp_filename)
            except OSError as err:
                print('*** Could not delete temp file "{0}" ***'.format(temp_filename))
                print(err)

    if file_data is None:
        # insufficient data for simulation
        return ERROR_TUPLE

    # compute pseudoperson tuple distributions prior to remapping data
    _compute_pseudoperson_distributions(file_data)

    # Remap the categorical varable values to consecutive integers.
    # Also do consistency checks on the remapped data.
    _remap_data(file_data)

    # compute empirical CDFs for each variable
    cdf_list = _compute_cdfs(variable_names, record_count)

    # compute Kendall's tau matrix for the remapped data
    data_list = _build_data_list(variable_names)
    tau = correlation_matrix.kendalls_tau_matrix(data_list)
    if tau is None:
        return ERROR_TUPLE

    return variable_names, tau, cdf_list, file_data


###############################################################################
def plot_marginal_distributions():
    """
    Generate plots of the marginal distributions for each variable.
    """

    plots_hl7.plot_marginals(_state)
    

###############################################################################
def plot_pdf_ecdf():
    """
    Generate plots of the set of ECDFs and the set of inverse ECDFs for
    each variable.
    """

    plots_hl7.plot_pdf_ecdf(_state)


###############################################################################
def get_remapped_data(variable_names):
    """
    A convenience function to return the remapped data arrays. This function
    is only needed for debugging and for checking matrix element convergence.
    """

    remapped_data = []
    
    for name in variable_names:
        if HL7.NAME_AGE == name:
            remapped_data.append(copy.deepcopy(_state.age_data))
        elif HL7.NAME_SEX == name:
            remapped_data.append(copy.deepcopy(_state.sex_data))
        elif HL7.NAME_RACE == name:
            remapped_data.append(copy.deepcopy(_state.race_data))
        elif HL7.NAME_ETHNICITY == name:
            remapped_data.append(copy.deepcopy(_state.ethnicity_data))
        elif HL7.NAME_COUNTY == name:
            remapped_data.append(copy.deepcopy(_state.county_data))
        elif HL7.NAME_CASE_STATUS == name:
            remapped_data.append(copy.deepcopy(_state.case_status_data))
        elif HL7.NAME_PREGNANT == name:
            remapped_data.append(copy.deepcopy(_state.pregnant_data))

    return remapped_data


###############################################################################
def remap_synthetic_results(X_list, variable_names):
    """
    Convert the synthetic data from remapped values to HL7 values.
    """

    assert len(X_list) == len(variable_names)
    
    syn_age = None
    syn_sex = None
    syn_race = None
    syn_ethnicity = None
    syn_county = None
    syn_case_status = None
    syn_pregnant = None
    
    for i,name in enumerate(variable_names):
        if HL7.NAME_AGE == name:
            syn_age = X_list[i]
        elif HL7.NAME_SEX == name:
            syn_sex = X_list[i]
        elif HL7.NAME_RACE == name:
            syn_race = X_list[i]
        elif HL7.NAME_ETHNICITY == name:
            syn_ethnicity = X_list[i]
        elif HL7.NAME_COUNTY == name:
            syn_county = X_list[i]
        elif HL7.NAME_CASE_STATUS == name:
            syn_case_status = X_list[i]
        elif HL7.NAME_PREGNANT == name:
            syn_pregnant = X_list[i]
            
    # remap all data
    remapped_age = []
    remapped_sex = []
    remapped_race = []
    remapped_ethnicity = []
    remapped_county = []
    remapped_case_status = []
    remapped_pregnant = []
    
    if syn_age is not None:
        for a in syn_age:
            if a in HL7.INV_AGE_MAP:
                remapped_age.append(HL7.INV_AGE_MAP[a])
            else:
                remapped_age.append(a)
    if syn_sex is not None:
        for s in syn_sex:
            if s in HL7.INV_SEX_MAP:
                remapped_sex.append(HL7.INV_SEX_MAP[s])
            else:
                remapped_sex.append(s)
    if syn_race is not None:
        for r in syn_race:
            if r in _state.inv_race_map:
                remapped_race.append(_state.inv_race_map[r])
            else:
                remapped_race.append(r)
    if syn_ethnicity is not None:
        for e in syn_ethnicity:
            if e in HL7.INV_ETHNICITY_MAP:
                remapped_ethnicity.append(HL7.INV_ETHNICITY_MAP[e])
            else:
                remapped_ethnicity.append(e)
    if syn_county is not None:
        for c in syn_county:
            if c in _state.inv_county_map:
                remapped_county.append(_state.inv_county_map[c])
            else:
                remapped_county.append(c)
    if syn_case_status is not None:
        for c in syn_case_status:
            if c in HL7.INV_CASE_STATUS_MAP:
                remapped_case_status.append(HL7.INV_CASE_STATUS_MAP[c])
            else:
                remapped_case_status.append(c)
    if syn_pregnant is not None:
        for p in syn_pregnant:
            if p in HL7.INV_PREGNANT_MAP:
                remapped_pregnant.append(HL7.INV_PREGNANT_MAP[p])
            else:
                remapped_pregnant.append(p)

    remapped_results = []
    for name in variable_names:
        if HL7.NAME_AGE == name:
            remapped_results.append(remapped_age)
        elif HL7.NAME_SEX == name:
            remapped_results.append(remapped_sex)
        elif HL7.NAME_RACE == name:
            remapped_results.append(remapped_race)
        elif HL7.NAME_ETHNICITY == name:
            remapped_results.append(remapped_ethnicity)
        elif HL7.NAME_COUNTY == name:
            remapped_results.append(remapped_county)
        elif HL7.NAME_CASE_STATUS == name:
            remapped_results.append(remapped_case_status)
        elif HL7.NAME_PREGNANT == name:
            remapped_results.append(remapped_pregnant)

    return remapped_results


###############################################################################
def is_valid_pseudoperson(sex_float, cs_float):

    pp_cat = _to_pseudoperson_category(sex_float, cs_float)
    return pp_cat in _state.pseudoperson_map


###############################################################################
def _select_date_tuple(sex_hl7,
                       case_status_hl7,
                       pp_categories,
                       pp_unexpected_map,
                       ref_dt,
                       max_dt_orig,
                       rng):
    """
    Given the HL7 values for sex and case status, compute the pseudoperson
    type and select a date tuple from the appropriate distribution.
    
    The array 'pp_categories' contains the keys of the pseudoperson map,
    which are the pseudoperson categories in the source data.
    """

    # max days allowed at the right edge of the dataset
    MAX_DAYS = 30

    # get the remapped floating point values for these
    sex_float = HL7.SEX_MAP[sex_hl7]
    cs_float  = HL7.CASE_STATUS_MAP[case_status_hl7]

    # find the pseudoperson category and get the object
    pp_cat = _to_pseudoperson_category(sex_float, cs_float)

    # All invalid pseudopersons should have been corrected prior to this
    # point. If pp_cat is a pseudoperson category that is NOT in the
    # pseudoperson map, just choose a random pseudoperson. This is not
    # ideal, but it should actually happen only rarely.
    if pp_cat not in _state.pseudoperson_map:
        # choose a random pseudoperson category
        #print('not in pseudoperson map: "{0}"'.format(pp_cat))
        pp_unexpected_map[pp_cat] += 1
        pp_cat = rng.choice(pp_categories)

    assert pp_cat in _state.pseudoperson_map
    pp_obj = _state.pseudoperson_map[pp_cat]

    # Select a date tuple whose max date from the reference date is
    # MAX_DAYS or fewer. Make repeated attempts; if all fail, truncate
    # the tuple at MAX_DAYS from max_dt_orig.

    success = False
    for i in range(1000):

        # generate a uniform random number on [0,1]
        u = rng.uniform(0.0, 1.0)

        # transform with the inverse ECDF for this tuple
        # distribution to get a synthetic tuple index
        tup_index = pp_obj.tuple_ecdf.inv(u)

        # lookup the tuple at this index
        date_tup = pp_obj.tuple_index_map[tup_index]

        # get max tuple component
        components = [c for c in date_tup if c is not None]
        if 0 == len(components):
            success = True
            break
        else:
            max_c = max(components)
            max_dt = ref_dt + datetime.timedelta(days=max_c)
            if max_dt <= max_dt_orig:
                success = True
                break

    if not success:
        # truncate entire tuple at max_dt_orig
        new_tuple = []
        for c in date_tup:
            if c is None or 0 == c:
                new_tuple.append(c)
            else:
                # truncate this component if beyond the date limit
                max_day_diff = (max_dt_orig - ref_dt).days + MAX_DAYS
                new_tuple.append(min(c, max_day_diff))
        date_tup = tuple(new_tuple)

    return date_tup
    

###############################################################################
def generate_date_tuples(num_samples,
                         signal,
                         dates,
                         variable_names,
                         synthetic_results,
                         str_max_date_orig,
                         rng):
    """
    Generate the date tuples for the synthetic data. One distinct tuple per
    synthetic sample is required.
    """

    # convert the max date to a datetime object, needed for date arithmetic
    max_dt_orig = datetime.datetime.strptime(str_max_date_orig, '%Y-%m-%d')

    signal_len = len(signal)

    # keeps track of synthetic pseudopersons not in the orignal data
    pp_unexpected_map = defaultdict(int)

    # pseudoperson categories
    pp_categories = [k for k in _state.pseudoperson_map.keys()]

    # synthetic sex and case status, needed to determine pseudoperson category
    sex_s = []
    case_status_s = []

    for i,name in enumerate(variable_names):
        if HL7.NAME_SEX == name:
            sex_s = synthetic_results[i]
        elif HL7.NAME_CASE_STATUS == name:
            case_status_s = synthetic_results[i]

    # variable 'i' indexes entries in 'signal' and 'dates'
    # range: [0, len(signal))
    i = 0

    # variable 'j' indexes the synthetic data points
    # range: [0, num_samples)
    j = 0

    # ref_date is in 'YYYY-mm-dd' format
    ref_date = dates[0]
    ref_datetime = datetime.datetime.strptime(ref_date, '%Y-%m-%d')

    # list of generated date tuples
    # 0-count dates receive an entry of None
    date_tuple_list = []

    # boolean indicating whether all possible dates have been generated
    exhausted_all_dates = False

    while True:
        # number of synthetic case reports on ref_date
        count = int(signal[i % signal_len])
        assert count >= 0

        if 0 == count:
            # no synthetic sample is consumed for 0-count entries, and no
            # date tuple needs to be generated; do not increment j
            pass
        else:
            # generate 'count' synthetic date tuples, all on this ref_date
            for k in range(count):
                sex         = _EMPTY_STRING
                case_status = _EMPTY_STRING
                
                # get sex and case_status for the next synthetic sample
                if len(sex_s) > 0:
                    sex = sex_s[j]
                if len(case_status_s) > 0:
                    case_status = case_status_s[j]

                # Generate a date tuple. Note that the sex and case_status
                # values have been converted to HL7 by this point.
                if _EMPTY_STRING != sex and _EMPTY_STRING != case_status:
                    assert sex in HL7.SEX_MAP
                    assert case_status in HL7.CASE_STATUS_MAP
                    tup = _select_date_tuple(sex,
                                             case_status,
                                             pp_categories,
                                             pp_unexpected_map,
                                             ref_datetime,
                                             max_dt_orig,
                                             rng)
                else:
                    # no date tuple; only a single field assigned the ref_date
                    tup = None

                date_tuple_list.append(tup)

                # that's one more tuple generated
                j += 1

                # finished if 'num_samples' have been generated
                if j >= num_samples:
                    break

        # go to the next date and count
        i += 1

        # finished if 'num_samples' have been written
        if j >= num_samples:
            break

        # keep going into the future if past the end of the signal
        if i >= signal_len:
            # 'cur_date' is the ref_datetime just written
            cur_date = ref_datetime

            # the next ref_date is one day later
            if cur_date.year == datetime.MAXYEAR:
                if 12 == cur_date.month and 31 == cur_date.day:
                    # exhausted all dates, so stop here
                    exhausted_all_dates = True
                    break

            next_day = cur_date + datetime.timedelta(days=1)
            ref_date = next_day.strftime('%Y-%m-%d')
        else:
            ref_date = dates[i]

        ref_datetime = datetime.datetime.strptime(ref_date, '%Y-%m-%d')        

    if len(pp_unexpected_map) > 0:
        print('Synthetic pseudopersons NOT in original data: ')
        for pp_cat, count in pp_unexpected_map.items():
            symbol = HL7.PSEUDOPERSON_SYMBOL_MAP[pp_cat]
            print('\t{0} => {1}'.format(symbol, count))
    
    # should have as many tuple entries as synthetic samples
    # (some entries could be None)
    print('Generated {0} date tuples for {1} synthetic samples.'.
          format(len(date_tuple_list), num_samples))
    if not exhausted_all_dates and len(date_tuple_list) != num_samples:
        raise SystemExit('generate_date_tuples')

    return date_tuple_list


###############################################################################
def _output_tup_to_string(output_tuple):
    """
    Generate a string to be written into the output file, csv format only.
    """
    
    # take the fields in order and concatenate
    fields = [str(f) for f in output_tuple]
    line = ','.join(fields)
    return line


###############################################################################
def _to_zero_tuple(ref_date):
    """
    Assign the reference date to the appropriate field in the output tuple.
    """

    if HL7.FIELD_REF_DATE == HL7.FIELD_DATE1:
        # first_elec_submit_dt is the reference date
        zero_tup = _OutputTuple(
            first_elec_submit_dt = ref_date,
            count = '0'
            # other fields default to the empty string
        )
    else:
        # report_dt_new is the reference date
        zero_tup = _OutputTuple(
            report_dt_new = ref_date,
            count = '0'
            # other fields default to the empty string
        )

    return zero_tup
    

###############################################################################
def _to_zero_dict(ref_date):
    """
    Assign the reference date to the appropriate field in the output record.
    """

    if HL7.FIELD_REF_DATE == HL7.FIELD_DATE1:
        # first_elec_submit_dt is the reference date
        zero_record = {
            HL7.FIELD_DATE1:ref_date,
            HL7.FIELD_COUNT:0
        }
    else:
        # report_dt_new is the reference date
        zero_record = {
            HL7.FIELD_REPORT_DATE:ref_date,
            HL7.FIELD_COUNT:0
        }

    return zero_record


###############################################################################
def write_csv_file(filename,
                   num_samples,
                   signal,
                   dates,
                   variable_names,
                   synthetic_results,
                   max_dt_orig,
                   date_tuple_list):
    """
    This function writes the synthetic data to a CSV file.
    """
    
    if _TRACE:
        print('From write_csv_file...')
        print('\tnum_samples: {0}'.format(num_samples))
        print('\tlen(synthetic_results[0]): {0}'.format(len(synthetic_results[0])))

    assert len(synthetic_results[0]) >= num_samples

    signal_len = len(signal)

    # synthetic data
    age_s = []
    sex_s = []
    race_s = []
    ethnicity_s = []
    case_status_s = []
    county_s = []
    pregnant_s = []

    for i,name in enumerate(variable_names):
        if HL7.NAME_AGE == name:
            age_s = synthetic_results[i]
        elif HL7.NAME_SEX == name:
            sex_s = synthetic_results[i]
        elif HL7.NAME_RACE == name:
            race_s = synthetic_results[i]
        elif HL7.NAME_ETHNICITY == name:
            ethnicity_s = synthetic_results[i]
        elif HL7.NAME_CASE_STATUS == name:
            case_status_s = synthetic_results[i]
        elif HL7.NAME_COUNTY == name:
            county_s = synthetic_results[i]
        elif HL7.NAME_PREGNANT == name:
            pregnant_s = synthetic_results[i]

    with open(filename, 'wt') as outfile:
        csv_text = _output_tup_to_string(tuple(HL7.OUTPUT_FIELDS))
        outfile.write('{0}\n'.format(csv_text))

        # variable 'i' indexes entries in 'signal' and 'dates'
        # range: [0, len(signal))
        i = 0

        # variable 'j' indexes the synthetic data points
        # range: [0, num_samples)
        j = 0

        # ref_date is in 'YYYY-mm-dd' format
        ref_date = dates[0]
        ref_datetime = datetime.datetime.strptime(ref_date, '%Y-%m-%d')
        while True:
            count = int(signal[i % signal_len])
            assert count >= 0

            if 0 == count:
                # special handling for zero count - do not increment j, since
                # no synthetic sample is consumed
                zero_tup = _to_zero_tuple(ref_date)
                csv_text = _output_tup_to_string(zero_tup)
                outfile.write('{0}\n'.format(csv_text))
            else:
                # insert 'count' synthetic tuples, same date, all have count=1
                for k in range(count):

                    age         = _EMPTY_STRING
                    age_units   = _EMPTY_STRING
                    sex         = _EMPTY_STRING
                    race        = _EMPTY_STRING
                    ethnicity   = _EMPTY_STRING
                    case_status = _EMPTY_STRING
                    county      = _EMPTY_STRING
                    pregnant    = _EMPTY_STRING

                    if len(age_s) > 0:
                        age = int(age_s[j])
                        age_units = HL7.AGE_UNITS_YEARS
                    if len(sex_s) > 0:
                        sex = sex_s[j]
                    if len(race_s) > 0:
                        race = race_s[j]
                    if len(ethnicity_s) > 0:
                        ethnicity = ethnicity_s[j]
                    if len(case_status_s) > 0:
                        case_status = case_status_s[j]
                    if len(county_s) > 0:
                        county = county_s[j]
                    if len(pregnant_s) > 0:
                        pregnant = pregnant_s[j]

                    # if male, set 'pregnant' field to empty string (not unk)
                    if HL7.SEX_MALE == sex:
                        pregnant = _EMPTY_STRING
                    
                    report_dt_new        = _EMPTY_STRING
                    first_elec_submit_dt = _EMPTY_STRING
                    diag_dt              = _EMPTY_STRING
                    died_dt              = _EMPTY_STRING
                    hosp_admit_dt        = _EMPTY_STRING
                    illness_onset_dt     = _EMPTY_STRING
                    invest_start_dt      = _EMPTY_STRING

                    # select the next date tuple, which encodes the relevant dates
                    tup = date_tuple_list[j]

                    # the ref_datetime matches the anchor date for all tuples
                    # (either the earliest date or the median date)

                    # compute datetimes for all correlated dates
                    if tup is not None and tup[0] is not None:
                        report_dt_new = ref_datetime + datetime.timedelta(days=tup[0])
                    if tup is not None and tup[1] is not None:
                        first_elec_submit_dt = ref_datetime + datetime.timedelta(days=tup[1])
                    if tup is not None and tup[2] is not None:
                        diag_dt = ref_datetime + datetime.timedelta(days=tup[2])
                    if tup is not None and tup[3] is not None:
                        died_dt = ref_datetime + datetime.timedelta(days=tup[3])
                    if tup is not None and tup[4] is not None:
                        hosp_admit_dt = ref_datetime + datetime.timedelta(days=tup[4])
                    if tup is not None and tup[5] is not None:
                        illness_onset_dt = ref_datetime + datetime.timedelta(days=tup[5])
                    if tup is not None and tup[6] is not None:
                        invest_start_dt = ref_datetime + datetime.timedelta(days=tup[6])

                    # assign the ref date to first_elec_submit_dt if no adequate tuple
                    if _EMPTY_STRING == first_elec_submit_dt:
                        first_elec_submit_dt = ref_datetime

                    output_tup = _OutputTuple(
                        report_dt_new         = report_dt_new,
                        count                 = '1',
                        age                   = age,
                        age_units             = age_units,
                        sex                   = sex,
                        ethnicity_txt         = ethnicity,
                        race_mapped           = race,
                        case_status_txt       = case_status,
                        notif_result_status   = HL7.NOTIF_RESULT_STATUS_FINAL,
                        pregnant              = pregnant,
                        first_elec_submit_dt  = first_elec_submit_dt,
                        subj_county           = county,
                        diag_dt               = diag_dt,
                        died_dt               = died_dt,
                        hosp_admit_dt         = hosp_admit_dt,
                        illness_onset_dt      = illness_onset_dt,
                        invest_start_dt       = invest_start_dt
                    )
                    csv_text = _output_tup_to_string(output_tup)
                    outfile.write('{0}\n'.format(csv_text))
                    j += 1

                    # finished if 'num_samples' have been written
                    if j >= num_samples:
                        break

            # go to the next date and count
            i += 1

            # finished if 'num_samples' have been written
            if j >= num_samples:
                break

            # keep going into the future if past the end of the signal
            if i >= signal_len:
                # 'cur_date' is the ref_datetime just written
                cur_date = ref_datetime

                # the next ref_date is one day later
                if cur_date.year == datetime.MAXYEAR:
                    if 12 == cur_date.month and 31 == cur_date.day:
                        # exhausted all dates
                        print('*** reached max date of 9999-12-31 ***')
                        print('output file contains {0} samples'.format(j))
                        return

                next_day = cur_date + datetime.timedelta(days=1)
                ref_date = next_day.strftime('%Y-%m-%d')
            else:
                ref_date = dates[i]

            ref_datetime = datetime.datetime.strptime(ref_date, '%Y-%m-%d')

        # append 0-count entries if not yet to the end of the signal
        if i < signal_len:
            while i < signal_len:
                ref_date = dates[i]
                zero_tup = _to_zero_tuple(ref_date)
                csv_text = _output_tup_to_string(zero_tup)
                outfile.write('{0}\n'.format(csv_text))
                i += 1


###############################################################################
def _terminate_json_file(outfile):

    # remove final comma and newline (both are single-byte chars in
    # ascii and utf-8); replace with end bracket
    current_pos = outfile.tell()
    outfile.seek(current_pos-2, os.SEEK_SET)
    outfile.truncate()
    outfile.write('\n]\n')


###############################################################################
def write_json_file(filename,
                    num_samples,
                    signal,
                    dates,
                    variable_names,
                    synthetic_results,
                    max_dt_orig,
                    date_tuple_list):
    """
    This function writes the synthetic data to a JSON file.
    """

    if _TRACE:
        print('From write_json_file...')
        print('\tnum_samples: {0}'.format(num_samples))
        print('\tlen(synthetic_results[0]): {0}'.format(len(synthetic_results[0])))

    assert len(synthetic_results[0]) >= num_samples

    signal_len = len(signal)

    # synthetic data
    age_s         = []
    sex_s         = []
    race_s        = []
    ethnicity_s   = []
    case_status_s = []
    county_s      = []
    pregnant_s    = []

    for i,name in enumerate(variable_names):
        if HL7.NAME_AGE == name:
            age_s = synthetic_results[i]
        elif HL7.NAME_SEX == name:
            sex_s = synthetic_results[i]
        elif HL7.NAME_RACE == name:
            race_s = synthetic_results[i]
        elif HL7.NAME_ETHNICITY == name:
            ethnicity_s = synthetic_results[i]
        elif HL7.NAME_CASE_STATUS == name:
            case_status_s = synthetic_results[i]
        elif HL7.NAME_COUNTY == name:
            county_s = synthetic_results[i]
        elif HL7.NAME_PREGNANT == name:
            pregnant_s = synthetic_results[i]

    with open(filename, 'wt') as outfile:
        outfile.write('[\n')

        # variable 'i' indexes entries in 'signal' and 'dates'
        # range: [0, len(signal))
        i = 0

        # variable 'j' indexes the synthetic data points
        # range: [0, num_samples)
        j = 0
        
        # ref_date is in 'YYYY-mm-dd' format
        ref_date = dates[0]
        ref_datetime = datetime.datetime.strptime(ref_date, '%Y-%m-%d')
        while True:
            count = int(signal[i % signal_len])
            assert count >= 0

            if 0 == count:
                # special handling for zero count - do not increment j, since
                # no synthetic sample is consumed
                record = _to_zero_dict(ref_date)
                json.dump(record, outfile, indent=4)
                outfile.write(',\n')
            else:
                # insert 'count' synthetic tuples, same date, all have count=1
                for k in range(count):

                    record = {
                        HL7.FIELD_COUNT:1
                    }

                    if len(age_s) > 0:
                        record[HL7.FIELD_AGE] = int(age_s[j])
                        record[HL7.FIELD_AGE_UNITS] = HL7.AGE_UNITS_YEARS
                    if len(sex_s) > 0:
                        record[HL7.FIELD_SEX] = sex_s[j]
                    if len(race_s) > 0:
                        record[HL7.FIELD_RACE] = race_s[j]
                    if len(ethnicity_s) > 0:
                        record[HL7.FIELD_ETHNICITY] = ethnicity_s[j]
                    if len(case_status_s) > 0:
                        record[HL7.FIELD_CASE_STATUS] = case_status_s[j]
                    if len(county_s) > 0:
                        record[HL7.FIELD_SUBJ_COUNTY] = county_s[j]
                    if len(pregnant_s) > 0:
                        record[HL7.FIELD_PREGNANT] = pregnant_s[j]

                    # if male, set 'pregnant' field to empty string (not unk)
                    if HL7.SEX_MALE == record[HL7.FIELD_SEX]:
                        record[HL7.FIELD_PREGNANT] = _EMPTY_STRING
                        
                    report_dt_new        = _EMPTY_STRING
                    first_elec_submit_dt = _EMPTY_STRING
                    diag_dt              = _EMPTY_STRING
                    died_dt              = _EMPTY_STRING
                    hosp_admit_dt        = _EMPTY_STRING
                    illness_onset_dt     = _EMPTY_STRING
                    invest_start_dt      = _EMPTY_STRING

                    # select the next date tuple, which encodes the relevant dates                    
                    tup = date_tuple_list[j]

                    # the ref_datetime matches the earliest date for all tuples
                    # (either the earliest date or the median date)

                    # compute datetimes for all correlated dates
                    if tup is not None and tup[0] is not None:
                        report_dt_new = ref_datetime + datetime.timedelta(days=tup[0])
                        report_dt_new = str(report_dt_new)
                    if tup is not None and tup[1] is not None:
                        first_elec_submit_dt = ref_datetime + datetime.timedelta(days=tup[1])
                        first_elec_submit_dt = str(first_elec_submit_dt)
                    if tup is not None and tup[2] is not None:
                        diag_dt = ref_datetime + datetime.timedelta(days=tup[2])
                        diag_dt = str(diag_dt)
                    if tup is not None and tup[3] is not None:
                        died_dt = ref_datetime + datetime.timedelta(days=tup[3])
                        died_dt = str(died_dt)
                    if tup is not None and tup[4] is not None:
                        hosp_admit_dt = ref_datetime + datetime.timedelta(days=tup[4])
                        hosp_admit_dt = str(hosp_admit_dt)
                    if tup is not None and tup[5] is not None:
                        illness_onset_dt = ref_datetime + datetime.timedelta(days=tup[5])
                        illness_onset_dt = str(illness_onset_dt)
                    if tup is not None and tup[6] is not None:
                        invest_start_dt = ref_datetime + datetime.timedelta(days=tup[6])
                        invest_start_dt = str(invest_start_dt)

                    # assign the ref date to first_elec_submit_dt if no adequate tuple
                    if _EMPTY_STRING == first_elec_submit_dt:
                        first_elec_submit_dt = str(ref_datetime)

                    record[HL7.FIELD_REPORT_DATE] = report_dt_new
                    record[HL7.FIELD_DATE1]       = first_elec_submit_dt
                    record[HL7.FIELD_DATE2]       = diag_dt
                    record[HL7.FIELD_DATE3]       = died_dt
                    record[HL7.FIELD_DATE6]       = hosp_admit_dt
                    record[HL7.FIELD_DATE7]       = illness_onset_dt
                    record[HL7.FIELD_DATE8]       = invest_start_dt
                    record[HL7.FIELD_NRS]         = HL7.NOTIF_RESULT_STATUS_FINAL

                    json.dump(record, outfile, indent=4)
                    outfile.write(',\n')
                    j += 1

                    # finished if 'num_samples' have been written
                    if j >= num_samples:
                        break

            # go to the next date and count
            i += 1

            # finished if 'num_samples' have been written
            if j >= num_samples:
                break

            # keep going into the future if past the end of the signal
            if i >= signal_len:
                # 'cur_date' is the ref_date just written
                cur_date = ref_datetime

                # the next ref_date is one day later
                if cur_date.year == datetime.MAXYEAR:
                    if 12 == cur_date.month and 31 == cur_date.day:
                        # exhausted all dates
                        print('*** reached max date of 9999-12-31 ***')
                        print('output file contains {0} samples'.format(j))
                        _terminate_json_file(outfile)
                        return

                next_day = cur_date + datetime.timedelta(days=1)
                ref_date = next_day.strftime('%Y-%m-%d')
            else:
                ref_date = dates[i]

            ref_datetime = datetime.datetime.strptime(ref_date, '%Y-%m-%d')

        # append 0-count entries if not yet to the end of the signal
        if i < signal_len:
            while i < signal_len:
                ref_date = dates[i]
                record = _to_zero_dict(ref_date)
                json.dump(record, outfile, indent=4)
                outfile.write(',\n')
                i += 1

        _terminate_json_file(outfile)


###############################################################################
def write_output_file(filename,
                      num_samples,
                      signal,
                      dates,
                      variable_names,
                      synthetic_results,
                      str_max_date_orig,
                      date_tuple_list):
    
    fullname, ext = os.path.splitext(filename)
    assert ext is not None
    ext = ext.lower()

    max_dt = datetime.datetime.strptime(str_max_date_orig, '%Y-%m-%d')

    if '.json' == ext:
        write_json_file(filename,
                        num_samples,
                        signal,
                        dates,
                        variable_names,
                        synthetic_results,
                        max_dt,
                        date_tuple_list)
    else:
        write_csv_file(filename,
                       num_samples,
                       signal,
                       dates,
                       variable_names,
                       synthetic_results,
                       max_dt,
                       date_tuple_list)
