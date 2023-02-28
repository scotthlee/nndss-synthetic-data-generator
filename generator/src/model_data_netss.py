import os
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

from . import plots
from . import netss as NETSS
from . import correlation_matrix

from .ecdf import EmpiricalCDF
from .kernel_density_estimation import KdeECDF

_MISSING_AGE       = 'missing_age_count'
_MISSING_AGETYPE   = 'missing_agetype_count'
_MISSING_SEX       = 'missing_sex_count'
_MISSING_RACE      = 'missing_race_count'
_MISSING_HISPANIC  = 'missing_hispanic_count'
_MISSING_COUNTY    = 'missing_county_count'
_MISSING_CASSTAT   = 'missing_casstat_count'
_INVALID_SEX       = 'invalid_sex_count'
_INVALID_RACE      = 'invalid_race_count'
_INVALID_HISPANIC  = 'invalid_hispanic_count'
_INVALID_CASSTAT   = 'invalid_casstat_count'
_UNHANDLED_AGETYPE = 'unhandled_agetype_count'

# set to True to display debug info
_TRACE = False

_NETSS_RECORD_FIELDS = [
    'event_date',
    'count',
    'age',
    'sex',
    'race',
    'hispanic',
    'casstat',
    'county'
]
_NetssRecord = namedtuple('_NetssRecord', _NETSS_RECORD_FIELDS)

# a class to hold all of the state variables
class State:

    def reset(self):
        # clear and initialize all state variables

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
        self.ecdf_race = None

        # ethnicity (hispanic)
        self.hispanic_data = []
        self.hispanic_dict = defaultdict(int)
        self.ecdf_hispanic = None

        # county
        self.county_data = []
        self.county_dict = defaultdict(int)
        self.county_map = {}
        self.inv_county_map = {}
        self.county_max = None
        self.ecdf_county = None

        # case status (casstat)
        self.casstat_data = []
        self.casstat_dict = defaultdict(int)
        self.ecdf_casstat = None

    def __init__(self):
        self.reset()

# instance of the state class
_state = State()


###############################################################################
def enable_debug():
    global _TRACE
    _TRACE = True


###############################################################################
def _merge_files(filepath_list, rng):
    """
    Merge the data in the filepath list and write the merged data to
    a temp file.

    All files in the list are assumed to exist. The front-end code checks each
    file for existence.

    The file data is assumed to contain a header line followed by one or more
    lines with the following fields:

        EVENTD,COUNT,AGE,AGETYPE,SEX,RACE,HISPANIC,CASSTAT,COUNTY

    The first entry in each line is the event date in YYYY-MM-DD format.
    """

    DATE_FORMAT = '%Y-%m-%d'
    HEADER = 'EVENTD,COUNT,AGE,AGETYPE,SEX,RACE,HISPANIC,CASSTAT,COUNTY'

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
                    continue
                text = line.strip()
                if 0 == len(text):
                    continue

                # get components from the next line of data
                str_eventd, str_count, str_age, str_agetype, str_sex, str_race, \
                    str_hispanic, str_casstat, str_county = text.split(',')

                # update the start_date and record the final date seen
                if file_start_date is None:
                    file_start_date = datetime.datetime.strptime(str_eventd,
                                                                 DATE_FORMAT)
                    if start_date is None or file_start_date < start_date:
                        start_date = file_start_date

                # record the final date seen, which will be used to determine
                # the end date for this file
                file_end_date = str_eventd

                # if COUNT > 0 keep the line, otherwise discard
                if str_count is not None and len(str_count) > 0:
                    count = float(str_count)                
                    if 0 == count:
                        continue
                    else:
                        data_dict[str_eventd].append(text)
                else:
                    # assume this line is corrupted
                    continue

            # update the end_date
            end_date_this_file = datetime.datetime.strptime(file_end_date,
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

    #print('From _merge_files(): ')
    #print('      Start Date: {0}'.format(start_date))
    #print('        End Date: {0}'.format(end_date))
    #print('Merged File Name: {0}'.format(outfile_name))

    # write the merged data to disk, using COUNT==0 records for missing data
    cur_date = start_date
    with open(outfile_name, 'wt') as outfile:
        # write the header
        outfile.write('{0}\n'.format(HEADER))
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
                    outfile.write('{0}\n'.format(r))
            else:
                # write an empty record for the current date
                outfile.write('{0},0,,,,,,,\n'.format(str_cur_date))

            # increment to the next calendar day (handles leap years)
            cur_date += datetime.timedelta(days=1)
            
    return outfile_name


###############################################################################
def _empty_record(str_eventd):
    """
    Create and return an empty _NetssRecord for the given date.
    """
    
    return _NetssRecord(str_eventd, count=0, age=None, sex=None,
                        race=None, hispanic=None, casstat=None, county=None)


###############################################################################
def _to_netss_record(count, 
                     str_eventd,
                     str_count,
                     str_age,
                     str_agetype,
                     str_sex,
                     str_race,
                     str_hispanic,
                     str_casstat,
                     str_county):
    """
    Decode the data from a single line of an input file.
    """

    age      = None
    sex      = None
    race     = None
    hispanic = None
    county   = None
    casstat  = None
    problem_counts = defaultdict(int)

    # file data is stored as floats, so need float() to convert

    bad_age = str_age is None or 0 == len(str_age)
    if bad_age:
        # if the age is missing, substitute _AGE_UNKNOWN
        # can ignore AGETYPE in this case
        problem_counts[_MISSING_AGE] += 1
        age = float(NETSS.AGE_UNKNOWN)
    else:
        age = float(str_age)

        # check the agetype
        bad_agetype = str_agetype is None or 0 == len(str_agetype)
        if bad_agetype:
            # agetype is missing, so cannot interpret the age field
            problem_counts[_MISSING_AGETYPE] += 1
            age = float(NETSS.AGE_UNKNOWN)
        else:
            agetype = int(float(str_agetype))

            # convert age to years according to agetype
            if 0 == agetype:
                # years, no conversion needed
                pass
            elif 1 == agetype:
                # months
                age = age / 12.0
            elif 2 == agetype:
                # weeks
                age = age / 52.0
            elif 3 == agetype:
                # days
                age = age / 365.0
            elif NETSS.AGETYPE_UNKNOWN == agetype:
                # agetype is unknown, so cannot interpret the age field
                age = float(NETSS.AGE_UNKNOWN)
            else:
                print('Found unhandled agetype "{0}"'.format(agetype))
                problem_counts[_UNHANDLED_AGETYPE] += 1
                # treat as unknown for now; this is incorrect
                age = float(NETSS.AGE_UNKNOWN)

    if str_sex is None or 0 == len(str_sex):
        problem_counts[_MISSING_SEX] += 1
        sex = float(NETSS.SEX_UNKNOWN)
    else:
        sex = float(str_sex)
        if int(sex) not in NETSS.SEX_MAP:
            # invalid value
            problem_counts[_INVALID_SEX] += 1
            sex = float(NETSS.SEX_UNKNOWN)

    if str_race is None or 0 == len(str_race):
        problem_counts[_MISSING_RACE] += 1
        race = float(NETSS.RACE_UNKNOWN)
    else:
        race = float(str_race)
        if int(race) not in NETSS.RACE_MAP:
            # invalid value
            problem_counts[_INVALID_RACE] += 1
            race = float(NETSS.RACE_UNKNOWN)

    if str_hispanic is None or 0 == len(str_hispanic):
        problem_counts[_MISSING_HISPANIC] += 1
        hispanic = float(NETSS.HISPANIC_UNKNOWN)
    else:
        hispanic = float(str_hispanic)
        if int(hispanic) not in NETSS.HISPANIC_MAP:
            # invalid value
            problem_counts[_INVALID_HISPANIC] += 1
            hispanic = float(NETSS.HISPANIC_UNKNOWN)

    if str_casstat is None or 0 == len(str_casstat):
        problem_counts[_MISSING_CASSTAT] += 1
        casstat = float(NETSS.CASSTAT_UNKNOWN)
    else:
        casstat = float(str_casstat)
        if int(casstat) not in NETSS.CASSTAT_MAP:
            # invalid value
            problem_counts[_INVALID_CASSTAT] += 1
            casstat = float(NETSS.CASSTAT_UNKNOWN)

    if str_county is None or 0 == len(str_county):
        problem_counts[_MISSING_COUNTY] += 1
        county = float(NETSS.COUNTY_UNKNOWN)
    else:
        county = float(str_county)

    assert str_eventd is not None
    assert count is not None
    assert age is not None
    assert sex is not None
    assert race is not None
    assert hispanic is not None
    assert casstat is not None
    assert county is not None

    record = _NetssRecord(
        event_date = str_eventd,
        count      = count,
        age        = age,
        sex        = sex,
        race       = race,
        hispanic   = hispanic,
        casstat    = casstat,
        county     = county
    )

    return record, problem_counts


###############################################################################
def _load_single_file(filepath):
    """
    Data is assumed to contain these columns:

        EVENTD,COUNT,AGE,AGETYPE,SEX,RACE,HISPANIC,CASSTAT,COUNTY

    All values are floats except for EVENTD.
    
    EVENTD column: date in YYYY-MM-DD format
    COUNT column: either 1 or 0; ignore if 0, meaning no reports for that day
    AGE column: interpretation depends on the AGETYPE column
    AGETYPE column:  0 = years, range 0 to 120 years
                     1 = months, range 0-11 months
                     2 = weeks, range 0-52 weeks
                     3 = days, range 0-28 days
                     4 = census coded age groups
                     9 = unknown
    SEX column:      1 = male
                     2 = female
                     9 = unknown
    RACE column:     1 = Native American / Alaska Native
                     2 = Asian / Pacific Islander
                     3 = Afro American
                     5 = White
                     8 = Other
                     9 = Unknown
    HISPANIC column: 1 = Hispanic
                     2 = Not Hispanic
                     9 = Unknown
    CASSTAT column:  1 = Confirmed
                     2 = Probable
                     3 = Suspect
                     9 = Unknown
    COUNTY column:   FIPS code for reporting county, nonzero; the state code
                     (first two digits) is omitted
                     999 = Unknown
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
    total_problem_counts[_MISSING_AGETYPE] = 0
    total_problem_counts[_MISSING_SEX] = 0
    total_problem_counts[_MISSING_RACE] = 0
    total_problem_counts[_MISSING_HISPANIC] = 0
    total_problem_counts[_MISSING_COUNTY] = 0
    total_problem_counts[_MISSING_CASSTAT] = 0
    total_problem_counts[_INVALID_SEX] = 0
    total_problem_counts[_INVALID_RACE] = 0
    total_problem_counts[_INVALID_HISPANIC] = 0
    total_problem_counts[_INVALID_CASSTAT] = 0
    total_problem_counts[_UNHANDLED_AGETYPE] = 0

    records = []

    is_json_file = False
    if filepath.lower().endswith('csv'):

        with open(filepath, 'rt') as infile:
            # load file in single gulp
            raw = infile.read()
            lines = raw.split('\n')

            # skip header line
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

    for obj in records:
        if is_json_file:
            # 'obj' is a dict
            str_eventd   = None
            str_count    = None
            str_age      = None
            str_agetype  = None
            str_sex      = None
            str_race     = None
            str_hispanic = None
            str_casstat  = None
            str_county   = None

            if 'eventd' in obj:
                str_eventd = str(obj['eventd'])
            if 'count' in obj:
                str_count = str(obj['count'])
            if 'age' in obj:
                str_age = str(obj['age'])
            if 'agetype' in obj:
                str_agetype = str(obj['agetype'])
            if 'sex' in obj:
                str_sex = str(obj['sex'])
            if 'race' in obj:
                str_race = str(obj['race'])
            if 'hispanic' in obj:
                str_hispanic = str(obj['hispanic'])
            if 'casstat' in obj:
                str_casstat = str(obj['casstat'])
            if 'county' in obj:
                str_county = str(obj['county'])

        else:
            # obj is a string of comma-separated values
            # read data as strings and convert to int;
            # keep date as string in YYYY-MM-DD format
            str_eventd, str_count, str_age, str_agetype, str_sex, str_race, \
                str_hispanic, str_casstat, str_county = obj.split(',')

        if str_count is not None and len(str_count) > 0:
            count = float(str_count)                
            if 0 == count:
                zero_count_fields += 1
                file_data.append(_empty_record(str_eventd))
                continue
        else:
            # assume this line is corrupted
            bad_record_count += 1
            file_data.append(_empty_record(str_eventd))
            continue

        netss_record, problem_counts = _to_netss_record(count, 
                                                        str_eventd,
                                                        str_count,
                                                        str_age,
                                                        str_agetype,
                                                        str_sex,
                                                        str_race,
                                                        str_hispanic,
                                                        str_casstat,
                                                        str_county)

        for k,v in problem_counts.items():
            total_problem_counts[k] += v

        # count the number of different county values
        county_set.add(netss_record.county)

        # update some stats
        lines_with_data += 1
        count_sum += int(netss_record.count)

        file_data.append(netss_record)


    # required width for line counts
    line_count = len(records)
    width = len(str(line_count))
            
    print('Information for file "{0}":'.format(filepath))
    print('\t                   Start Date: {0}'.format(file_data[0].event_date))
    print('\t                     End Date: {0}'.format(file_data[-1].event_date))
    print('\tLine count (excluding header): {0}'.format(line_count, width))
    print('\t  Lines with complete records: {0:{1}d}'.format(lines_with_data, width))
    print('\t      Sum of all count fields: {0:{1}d}'.format(count_sum, width))
    print('\t        Lines with COUNT == 0: {0:{1}d}'.format(zero_count_fields, width))
    print('\t    No. of different counties: {0:{1}d}'.format(len(county_set), width))
    print('\t                  Bad records: {0:{1}d}'.format(bad_record_count, width))
    print('\t          Missing AGE entries: {0:{1}d}'.format(total_problem_counts[_MISSING_AGE],
                                                             width))
    print('\t      Missing AGETYPE entries: {0:{1}d}'.format(total_problem_counts[_MISSING_AGETYPE],
                                                             width))
    print('\t    Unhandled AGETYPE entries: {0:{1}d}'.format(total_problem_counts[_UNHANDLED_AGETYPE],
                                                             width))
    print('\t          Missing SEX entries: {0:{1}d}'.format(total_problem_counts[_MISSING_SEX],
                                                             width))
    print('\t         Missing RACE entries: {0:{1}d}'.format(total_problem_counts[_MISSING_RACE],
                                                             width))
    print('\t     Missing HISPANIC entries: {0:{1}d}'.format(total_problem_counts[_MISSING_HISPANIC],
                                                             width))
    print('\t       Missing COUNTY entries: {0:{1}d}'.format(total_problem_counts[_MISSING_COUNTY],
                                                             width))
    print('\t      Missing CASSTAT entries: {0:{1}d}'.format(total_problem_counts[_MISSING_CASSTAT],
                                                             width))
    print('\t          Invalid SEX entries: {0:{1}d}'.format(total_problem_counts[_INVALID_SEX],
                                                             width))
    print('\t         Invalid RACE entries: {0:{1}d}'.format(total_problem_counts[_INVALID_RACE],
                                                             width))
    print('\t     Invalid HISPANIC entries: {0:{1}d}'.format(total_problem_counts[_INVALID_HISPANIC],
                                                             width))
    print('\t      Invalid CASSTAT entries: {0:{1}d}'.format(total_problem_counts[_INVALID_CASSTAT],
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
def _get_counts(file_data):
    """
    Accumulate counts for each date, then return a list of (date, count) tuples.
    """
    
    date_count_dict = defaultdict(int)
    for r in file_data:
        date_count_dict[r.event_date] += r.count
        
    keys = sorted(date_count_dict.keys())
    tuple_list = [(k, date_count_dict[k]) for k in keys]
    return tuple_list


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
            hist_dict[i] = 0
            

###############################################################################
def _do_consistency_checks():
    """
    Check the data arrays to ensure that all arrays have the expected number
    of elements. Also check the PDFs.
    """

    age_counts = 0
    sex_counts = 0
    race_counts = 0
    hispanic_counts = 0
    county_counts = 0
    casstat_counts = 0

    for k in sorted(_state.age_dict.keys()):
        age_counts += _state.age_dict[k]
    for k in sorted(_state.sex_dict.keys()):
        sex_counts += _state.sex_dict[k]
    for k in sorted(_state.race_dict.keys()):
        race_counts += _state.race_dict[k]
    for k in sorted(_state.hispanic_dict.keys()):
        hispanic_counts += _state.hispanic_dict[k]
    for k in sorted(_state.county_dict.keys()):
        county_counts += _state.county_dict[k]    
    for k in sorted(_state.casstat_dict.keys()):
        casstat_counts += _state.casstat_dict[k]

    if _TRACE:
        print('AGE')
        for k in sorted(_state.age_dict.keys()):
            print('{0} => {1}'.format(k, _state.age_dict[k]))
        print('SEX')
        for k in sorted(_state.sex_dict.keys()):
            print('{0} => {1}'.format(k, _state.sex_dict[k]))
        print('RACE')
        for k in sorted(_state.race_dict.keys()):
            print('{0} => {1}'.format(k, _state.race_dict[k]))
        print('HISPANIC')
        for k in sorted(_state.hispanic_dict.keys()):
            print('{0} => {1}'.format(k, _state.hispanic_dict[k]))
        print('COUNTY')
        for k in sorted(_state.county_dict.keys()):
            print('{0} => {1}'.format(k, _state.county_dict[k]))
        print('CASSTAT')
        for k in sorted(_state.casstat_dict.keys()):
            print('{0} => {1}'.format(k, _state.casstat_dict[k]))

        print('     Sum of all sex counts: {0}'.format(sex_counts))
        print('    Sum of all race counts: {0}'.format(race_counts))
        print('Sum of all hispanic counts: {0}'.format(hispanic_counts))
        print('  Sum of all county counts: {0}'.format(county_counts))
        print(' Sum of all casstat counts: {0}'.format(casstat_counts))
        print('     Nonempty record count: {0}'.format(_state.nonempty_record_count))
        print('               _county_max: {0}'.format(_state.county_max))

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
    if hispanic_counts > 0:
        hispanic_pdf = _to_pdf(list(_state.hispanic_dict.values()))
        assert _state.nonempty_record_count == hispanic_counts
    if county_counts > 0:
        county_pdf = _to_pdf(list(_state.county_dict.values()))
        assert _state.nonempty_record_count == county_counts
    if casstat_counts > 0:
        casstat_pdf = _to_pdf(list(_state.casstat_dict.values()))
        assert _state.nonempty_record_count == casstat_counts

    if _TRACE:
        # check to see that the PDF bin values all sum to 1
        # (bins are contiguous with width == 1)
        sum_age_pdf = np.sum(age_pdf)
        print('     age_pdf bin sum: {0:.3f}'.format(sum_age_pdf))
        sum_sex_pdf = np.sum(sex_pdf)
        print('     sex_pdf bin sum: {0:.3f}'.format(sum_sex_pdf))
        sum_race_pdf = np.sum(race_pdf)
        print('    race_pdf bin sum: {0:.3f}'.format(sum_race_pdf))
        sum_hispanic_pdf = np.sum(hispanic_pdf)
        print('hispanic_pdf_bin_sum: {0:.3f}'.format(sum_hispanic_pdf))
        sum_county_pdf = np.sum(county_pdf)
        print('  county_pdf_bin_sum: {0:.3f}'.format(sum_county_pdf))
        sum_casstat_pdf = np.sum(casstat_pdf)
        print(' casstat_pdf_bin_sum: {0:.3f}'.format(sum_casstat_pdf))

    
###############################################################################
def _remap_data(file_data):
    """
    """

    # number of data entries with COUNT > 0
    _state.nonempty_record_count = 0

    for r in file_data:

        if r.count is None or 0 == r.count:
            continue

        # the file loader ensures that any line with a nonzero count has data
        _state.nonempty_record_count += 1

        # remap the age unknown value of 999 to -1
        if r.age is not None:
            if NETSS.AGE_UNKNOWN == r.age:
                the_age = NETSS.AGE_UNKNOWN_REMAP
            elif r.age > 0.0 and r.age < 1.0:
                # put fractions in the 0 bin
                the_age = 0.0
            elif r.age > 120.0:
                # 120 is the max age in years; put any values greater than this in the unknown bin
                the_age = NETSS.AGE_UNKNOWN_REMAP
            else:
                the_age = r.age

            _state.age_data.append(the_age)
            _state.age_dict[the_age] += 1

        # most of these variables have the unknown value remapped

        if r.sex is not None:
            assert r.sex in NETSS.SEX_MAP
            the_sex = NETSS.SEX_MAP[r.sex]
            _state.sex_data.append(the_sex)
            _state.sex_dict[the_sex] += 1

        if r.race is not None:
            assert r.race in NETSS.RACE_MAP
            the_race = NETSS.RACE_MAP[r.race]
            _state.race_data.append(the_race)
            _state.race_dict[the_race] += 1

        if r.hispanic is not None:
            assert r.hispanic in NETSS.HISPANIC_MAP
            the_hispanic = NETSS.HISPANIC_MAP[r.hispanic]
            _state.hispanic_data.append(the_hispanic)
            _state.hispanic_dict[the_hispanic] += 1

        if r.casstat is not None:
            assert r.casstat in NETSS.CASSTAT_MAP
            the_casstat = NETSS.CASSTAT_MAP[r.casstat]
            _state.casstat_data.append(the_casstat)
            _state.casstat_dict[the_casstat] += 1

        if r.county is not None:
            if NETSS.COUNTY_UNKNOWN == r.county:
                the_county = NETSS.COUNTY_UNKNOWN_REMAP
            else:
                the_county = r.county
            _state.county_data.append(the_county)
            _state.county_dict[the_county] += 1

    _zero_fill_missing_bins(NETSS.AGE_MIN, NETSS.AGE_MAX, _state.age_dict)
    _zero_fill_missing_bins(NETSS.SEX_MIN, NETSS.SEX_MAX, _state.sex_dict)
    _zero_fill_missing_bins(NETSS.RACE_MIN, NETSS.RACE_MAX, _state.race_dict)
    _zero_fill_missing_bins(NETSS.HISPANIC_MIN, NETSS.HISPANIC_MAX, _state.hispanic_dict)
    _zero_fill_missing_bins(NETSS.CASSTAT_MIN, NETSS.CASSTAT_MAX, _state.casstat_dict)
    # county handled separately below

    # remap the county variables to sequential integers
    if len(_state.county_data) > 0:
        new_county_dict = {}

        # take the keys in sorted order, 0 is now the unknown value
        keys = sorted(_state.county_dict.keys())

        index = 0
        for county_code in keys:
            # assign next county code to the next index
            _state.county_map[county_code] = index
            _state.inv_county_map[index] = county_code
            new_county_dict[index] = _state.county_dict[county_code]
            index += 1

        _state.county_dict = new_county_dict
        _state.county_max = len(_state.county_dict)

        # remap the county data to use the new values
        new_county_data = []
        for orig_county_code in _state.county_data:
            new_county_code = _state.county_map[orig_county_code]
            new_county_data.append(new_county_code)
        assert len(new_county_data) == len(_state.county_data)
        _state.county_data = new_county_data

    _do_consistency_checks()


###############################################################################
def _compute_cdfs(variable_names, record_count):
    """
    Compute the cumulative distribution functions for the data.
    """

    # use kernel density estimation for the age variable for smaller datasets
    if record_count < 2000:
        _state.ecdf_age = KdeECDF(_state.age_data,
                                  l=NETSS.AGE_MIN, r=NETSS.AGE_MAX, kernel='hybrid')
    else:
        _state.ecdf_age = EmpiricalCDF(_state.age_data)

    _state.ecdf_sex = EmpiricalCDF(_state.sex_data)
    _state.ecdf_race = EmpiricalCDF(_state.race_data)
    _state.ecdf_hispanic = EmpiricalCDF(_state.hispanic_data)
    _state.ecdf_county = EmpiricalCDF(_state.county_data)
    _state.ecdf_casstat = EmpiricalCDF(_state.casstat_data)

    ecdf_list = []

    # the order of the ECDF objects MUST match the order of the variable_names
    for name in variable_names:
        if NETSS.NAME_AGE == name:
            ecdf_list.append(_state.ecdf_age)
        elif NETSS.NAME_SEX == name:
            ecdf_list.append(_state.ecdf_sex)
        elif NETSS.NAME_RACE == name:
            ecdf_list.append(_state.ecdf_race)
        elif NETSS.NAME_HISPANIC == name:
            ecdf_list.append(_state.ecdf_hispanic)
        elif NETSS.NAME_COUNTY == name:
            ecdf_list.append(_state.ecdf_county)
        elif NETSS.NAME_CASSTAT == name:
            ecdf_list.append(_state.ecdf_casstat)

    return ecdf_list


###############################################################################
def _validate_names(variable_names):

    for name in variable_names:
        if name not in NETSS.NAME_STRING_SET:
            print('*** ERROR ***: invalid variable name "{0)"'.format(name))
            return False

    return True
        

###############################################################################
def _build_data_list(variable_names):
    """
    """

    data_list = []
    
    for name in variable_names:
        if NETSS.NAME_AGE == name:
            data_list.append(_state.age_data)
        elif NETSS.NAME_SEX == name:
            data_list.append(_state.sex_data)
        elif NETSS.NAME_RACE == name:
            data_list.append(_state.race_data)
        elif NETSS.NAME_HISPANIC == name:
            data_list.append(_state.hispanic_data)
        elif NETSS.NAME_COUNTY == name:
            data_list.append(_state.county_data)
        elif NETSS.NAME_CASSTAT == name:
            data_list.append(_state.casstat_data)

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
        try:
            os.remove(temp_filename)
        except OSError as err:
            print('*** Could not delete temp file "{0}" ***'.format(temp_filename))
            print(err)

    if file_data is None:
        # insufficient data for simulation
        return ERROR_TUPLE

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

    plots.plot_marginals(_state)


###############################################################################
def plot_pdf_ecdf():
    """
    Generate plots of the set of ECDFs and the set of inverse ECDFs for
    each variable.
    """

    plots.plot_pdf_ecdf(_state)
    

###############################################################################
def get_remapped_data(variable_names):
    """
    A convenience function to return the remapped data arrays. This function
    is only needed for debugging and for checking matrix element convergence.
    """

    remapped_data = []
    
    for name in variable_names:
        if NETSS.NAME_AGE == name:
            remapped_data.append(copy.deepcopy(_state.age_data))
        elif NETSS.NAME_SEX == name:
            remapped_data.append(copy.deepcopy(_state.sex_data))
        elif NETSS.NAME_RACE == name:
            remapped_data.append(copy.deepcopy(_state.race_data))
        elif NETSS.NAME_HISPANIC == name:
            remapped_data.append(copy.deepcopy(_state.hispanic_data))
        elif NETSS.NAME_COUNTY == name:
            remapped_data.append(copy.deepcopy(_state.county_data))
        elif NETSS.NAME_CASSTAT == name:
            remapped_data.append(copy.deepcopy(_state.casstat_data))

    return remapped_data
    

###############################################################################
def remap_synthetic_results(X_list, variable_names):
    """
    Convert the synthetic data from remapped values to NETSS values.
    """

    assert len(X_list) == len(variable_names)
    
    syn_age = None
    syn_sex = None
    syn_race = None
    syn_hispanic = None
    syn_county = None
    syn_casstat = None
    
    for i,name in enumerate(variable_names):
        if NETSS.NAME_AGE == name:
            syn_age = X_list[i]
        elif NETSS.NAME_SEX == name:
            syn_sex = X_list[i]
        elif NETSS.NAME_RACE == name:
            syn_race = X_list[i]
        elif NETSS.NAME_HISPANIC == name:
            syn_hispanic = X_list[i]
        elif NETSS.NAME_COUNTY == name:
            syn_county = X_list[i]
        elif NETSS.NAME_CASSTAT == name:
            syn_casstat = X_list[i]
            
    # remap all data
    remapped_age = []
    remapped_sex = []
    remapped_race = []
    remapped_hispanic = []
    remapped_county = []
    remapped_casstat = []
    
    if syn_age is not None:
        for a in syn_age:
            if a in NETSS.INV_AGE_MAP:
                remapped_age.append( NETSS.INV_AGE_MAP[a])
            else:
                remapped_age.append(a)
    if syn_sex is not None:
        for s in syn_sex:
            if s in NETSS.INV_SEX_MAP:
                remapped_sex.append( NETSS.INV_SEX_MAP[s])
            else:
                remapped_sex.append(s)
    if syn_race is not None:
        for r in syn_race:
            if r in NETSS.INV_RACE_MAP:
                remapped_race.append( NETSS.INV_RACE_MAP[r])
            else:
                remapped_race.append(r)
    if syn_hispanic is not None:
        for h in syn_hispanic:
            if h in NETSS.INV_HISPANIC_MAP:
                remapped_hispanic.append( NETSS.INV_HISPANIC_MAP[h])
            else:
                remapped_hispanic.append(h)
    if syn_county is not None:
        for c in syn_county:
            if c in _state.inv_county_map:
                remapped_county.append( _state.inv_county_map[c])
            else:
                remapped_county.append(c)
    if syn_casstat is not None:
        for c in syn_casstat:
            if c in NETSS.INV_CASSTAT_MAP:
                remapped_casstat.append( NETSS.INV_CASSTAT_MAP[c])
            else:
                remapped_casstat.append(c)

    remapped_results = []
    for name in variable_names:
        if NETSS.NAME_AGE == name:
            remapped_results.append(remapped_age)
        elif NETSS.NAME_SEX == name:
            remapped_results.append(remapped_sex)
        elif NETSS.NAME_RACE == name:
            remapped_results.append(remapped_race)
        elif NETSS.NAME_HISPANIC == name:
            remapped_results.append(remapped_hispanic)
        elif NETSS.NAME_COUNTY == name:
            remapped_results.append(remapped_county)
        elif NETSS.NAME_CASSTAT == name:
            remapped_results.append(remapped_casstat)

    return remapped_results


###############################################################################
def write_csv_file(filename,
                   num_samples,
                   signal,
                   dates,
                   variable_names,
                   synthetic_results):
    """
    This function writes the synthetic data to a CSV file.
    """

    HEADER = 'EVENTD,COUNT,AGE,AGETYPE,SEX,RACE,HISPANIC,CASSTAT,COUNTY'

    if _TRACE:
        print('From write_csv_file...')
        print('num_samples: {0}'.format(num_samples))
        print('len(synthetic_results[0]): {0}'.format(len(synthetic_results[0])))

    assert len(synthetic_results[0]) >= num_samples

    signal_len = len(signal)
    str_line = '{0},{1},{2},{3},{4},{5},{6},{7},{8}'

    # synthetic data
    age_s = []
    sex_s = []
    race_s = []
    hispanic_s = []
    casstat_s = []
    county_s = []

    for i,name in enumerate(variable_names):
        if NETSS.NAME_AGE == name:
            age_s = synthetic_results[i]
        elif NETSS.NAME_SEX == name:
            sex_s = synthetic_results[i]
        elif NETSS.NAME_RACE == name:
            race_s = synthetic_results[i]
        elif NETSS.NAME_HISPANIC == name:
            hispanic_s = synthetic_results[i]
        elif NETSS.NAME_CASSTAT == name:
            casstat_s = synthetic_results[i]
        elif NETSS.NAME_COUNTY == name:
            county_s = synthetic_results[i]

    with open(filename, 'wt') as outfile:
        outfile.write('{0}\n'.format(HEADER))

        # variable 'i' indexes entries in 'signal' and 'dates'
        # range: [0, len(signal))
        i = 0

        # variable 'j' indexes the synthetic data points
        # range: [0, num_samples)
        j = 0

        event_date = dates[0]
        while True:
            count = int(signal[i % signal_len])
            assert count >= 0

            if 0 == count:
                # special handling for zero count - do not increment j, since
                # no synthetic sample is consumed
                csv_text = str_line.format(event_date,0,'','','','','','','')
                outfile.write('{0}\n'.format(csv_text))
            else:
                # insert 'count' synthetic tuples, same date, all have count=1
                for k in range(count):

                    age = ''
                    agetype = ''
                    sex = ''
                    race = ''
                    hispanic = ''
                    casstat = ''
                    county = ''

                    if len(age_s) > 0:
                        age = int(age_s[j])
                        agetype = NETSS.AGETYPE_YEARS
                    if len(sex_s) > 0:
                        sex = int(sex_s[j])
                    if len(race_s) > 0:
                        race = int(race_s[j])
                    if len(hispanic_s) > 0:
                        hispanic = int(hispanic_s[j])
                    if len(casstat_s) > 0:
                        casstat = int(casstat_s[j])
                    if len(county_s) > 0:
                        county = int(county_s[j])

                    csv_text = str_line.format(event_date,
                                               1,
                                               age,
                                               agetype,
                                               sex,
                                               race,
                                               hispanic,
                                               casstat,
                                               county)
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
                # 'cur_date' is the event_date just written
                cur_date = datetime.datetime.strptime(event_date, '%Y-%m-%d')

                # the next event_date is one day later
                if cur_date.year == datetime.MAXYEAR:
                    if 12 == cur_date.month and 31 == cur_date.day:
                        # exhausted all dates
                        print('*** reached max date of 9999-12-31 ***')
                        print('output file contains {0} samples'.format(j))
                        return
                        
                next_day = cur_date + datetime.timedelta(days=1)
                event_date = next_day.strftime('%Y-%m-%d')
            else:
                event_date = dates[i]

        # append 0-count entries if not yet to the end of the signal
        if i < signal_len:
            while i < signal_len:
                event_date = dates[i]
                csv_text = str_line.format(event_date,0,'','','','','','','')
                outfile.write('{0}\n'.format(csv_text))
                i += 1


###############################################################################
def _terminate_json_file(outfile):
    
    # remove final comma and newline (both are single-byte chars in
    # ascii and utf-8); replace with end bracket
    current_pos = outfile.tell()
    outfile.seek(current_pos-2, os.SEEK_SET)
    outfile.truncate()
    outfile.write('\n]')
    

###############################################################################
def write_json_file(filename,
                    num_samples,
                    signal,
                    dates,
                    variable_names,
                    synthetic_results):
    """
    This function writes the synthetic data to a JSON file.
    """

    if _TRACE:
        print('From write_json_file...')
        print('num_samples: {0}'.format(num_samples))
        print('len(synthetic_results[0]): {0}'.format(len(synthetic_results[0])))

    assert len(synthetic_results[0]) >= num_samples

    signal_len = len(signal)

    # synthetic data
    age_s = []
    sex_s = []
    race_s = []
    hispanic_s = []
    casstat_s = []
    county_s = []

    for i,name in enumerate(variable_names):
        if NETSS.NAME_AGE == name:
            age_s = synthetic_results[i]
        elif NETSS.NAME_SEX == name:
            sex_s = synthetic_results[i]
        elif NETSS.NAME_RACE == name:
            race_s = synthetic_results[i]
        elif NETSS.NAME_HISPANIC == name:
            hispanic_s = synthetic_results[i]
        elif NETSS.NAME_CASSTAT == name:
            casstat_s = synthetic_results[i]
        elif NETSS.NAME_COUNTY == name:
            county_s = synthetic_results[i]

    # these two fields are always present
    KEY_DATE  = 'eventd'
    KEY_COUNT = 'count'

    with open(filename, 'wt') as outfile:
        outfile.write('[\n')

        # variable 'i' indexes entries in 'signal' and 'dates'
        # range: [0, len(signal))
        i = 0

        # variable 'j' indexes the synthetic data points
        # range: [0, num_samples)
        j = 0

        event_date = dates[0]
        while True:
            count = int(signal[i % signal_len])
            assert count >= 0

            if 0 == count:
                # special handling for zero count - do not increment j, since
                # no synthetic sample is consumed
                record = {
                    KEY_DATE:event_date,
                    KEY_COUNT:0
                }
                json.dump(record, outfile, indent=4)
                outfile.write(',\n')
            else:
                # insert 'count' synthetic tuples, same date, all have count=1
                for k in range(count):

                    record = {
                        KEY_DATE:event_date,
                        KEY_COUNT:1
                    }

                    if len(age_s) > 0:
                        record['age'] = int(age_s[j])
                        record['agetype'] = NETSS.AGETYPE_YEARS
                    if len(sex_s) > 0:
                        record['sex'] = int(sex_s[j])
                    if len(race_s) > 0:
                        record['race'] = int(race_s[j])
                    if len(hispanic_s) > 0:
                        record['hispanic'] = int(hispanic_s[j])
                    if len(casstat_s) > 0:
                        record['casstat'] = int(casstat_s[j])
                    if len(county_s) > 0:
                        record['county'] = int(county_s[j])

                    json.dump(record, outfile, indent=4)
                    outfile.write(',\n')

                    # finished if 'num_samples' have been written
                    j += 1
                    if j >= num_samples:
                        break

            # go to the next date and count
            i += 1

            # finished if 'num_samples' have been written
            if j >= num_samples:
                break

            # keep going into the future if past the end of the signal
            if i >= signal_len:
                # 'cur_date' is the event_date just written
                cur_date = datetime.datetime.strptime(event_date, '%Y-%m-%d')

                # the next event_date is one day later
                if cur_date.year == datetime.MAXYEAR:
                    if 12 == cur_date.month and 31 == cur_date.day:
                        # exhausted all dates
                        print('*** reached max date of 9999-12-31 ***')
                        print('output file contains {0} samples'.format(j))
                        _terminate_json_file(outfile)
                        return
                        
                next_day = cur_date + datetime.timedelta(days=1)
                event_date = next_day.strftime('%Y-%m-%d')
            else:
                event_date = dates[i]

        # append 0-count entries if not yet to the end of the signal
        if i < signal_len:
            while i < signal_len:
                event_date = dates[i]
                record = {
                    KEY_DATE:event_date,
                    KEY_COUNT:0
                }
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
                      synthetic_results):
    
    fullname, ext = os.path.splitext(filename)
    assert ext is not None
    ext = ext.lower()

    if '.json' == ext:
        write_json_file(filename,
                        num_samples,
                        signal,
                        dates,
                        variable_names,
                        synthetic_results)
    else:
        write_csv_file(filename,
                       num_samples,
                       signal,
                       dates,
                       variable_names,
                       synthetic_results)
