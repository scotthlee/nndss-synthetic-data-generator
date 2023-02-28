#!/usr/bin/env python
# coding: utf-8

# # create_preprocessed_data.ipynb
# This notebook contains code to lightly preprocess and organize the NETSS data (../gtarc.csv) and HL7 data (../MVPS/gen_core.csv) into folders such that the folder names represent jurisdictions and files within represent disease state codes. For example, NETSS data for varicella in Florida would be located under ./NETSS/Florida/10030.csv. The purpose is to make the synthetic data generation more efficient; instead of loading and processing the entire raw dataset each time, one can leverage the folder/file structure created from this script to load a ready-made subset directly.
# 
# The preprocessing performed is as follows:
# 
# NETSS Data:
# 1. Discard records where COUNT is not equal to 1, per CDC recommendations.
# 2. Discard records having a missing EVENTD.
# 3. Only retain records with complete years: Currently, 2014-2019, but generally 2014 through the year that this script is being run, minus 1.
# 4. For each STATE and EVENT type, create a zero-padded time series such that the series begins on the first day of the first month of the first event and ends on the last day of the last month of the last event.
# 
# HL7 Data (Gen Core):
# 1. Define a constant COUNT variable equal to 1 (to align with the NETSS 1 row = 1 COUNT methodology). This implies that 1 row = 1 record in the HL7 data.
# 2. Per a CSELS/GTRI meeting on 10/19/2020, subset to `current_record_flag == 'Y'` -- this is because trying to generate synthetic data consistently across multiple records is an order of magnitude more complicated and isn't perceived to yield much additional benefit.
# 3. Discard records having a missing report_date.
# 4. Only retain dates that have at least 1,000 total counts in the raw dataset (here, 2004-2020) -- others not meeting this criterion are considered anomalous (e.g. years 1951, 2109) or not having enough data to reliably make synthetic data (e.g. 2002, 2003).
# 5. For each `report_state` and `condition_code` type, create a zero-padded time series such that the series begins on the first day of the first month of the first event and ends on the last day of the last month of the last event.

# ## Setup

# In[1]:


# Imports
import collections
import datetime
import numpy as np
import os
import pandas as pd
import shutil
import sys
import time
import csv
import re


# In[2]:


# Start timer
start_time = time.time()


# In[3]:


# Define which columns will be retained from each dataset out of the possibilities
usecols_NETSS = [
 'AGE',
 'AGETYPE',
 #'BIRTHD',
 #'CASEID',
 'CASSTAT',
 #'CDCDATE',
 'COUNT',
 'COUNTY',
 #'COUNTY_OF_RESIDENCE',
 #'DATET',
 #'DISEASE',
 'EVENT',
 'EVENTD',
 #'EVENTN',
 #'EXPANDED_CASEID',
 #'EXT_DATE',
 #'GRAND_TOTAL',
 'HISPANIC',
 #'IMPORT',
 #'INT_DATE',
 #'INV501_cd',
 #'INV501_desc',
 #'NON_RES_COUNT',
 #'OTHER',
 #'OUTBR',
 #'PRINTED',
 'RACE',
 #'RACECAT',
 #'RECTYPE',
 'SEX',
 #'SITE',
 'STATE',
 #'STNAME',
 #'TER_COUNT',
 #'UPDATE',
 #'USA_COUNT',
 #'WEEK',
 #'WSYSTEM',
 #'YEAR',
 #'message_received',
 #'msg_seq_id',
 #'msg_transaction_id'
]


# In[4]:


usecols_HL7 = [
'age',
'age_units',
#'Age_Units_txt',
#'Birth_Country',
#'Birth_Country_oth',
#'Birth_Date',
'birth_date_str',
#'Case_Status',
'case_status_txt',
#'Cntry_Usual_Res',
#'Cntry_Usual_Res_oth',
#'Cntry_Usual_Res_txt',
#'Comment',
'condition_code',
#'Condition_Code_txt',
#'Condition_MMG',
#'Core_Uid',
'diag_dt',
#'Diag_Dt_str',
#'Died',
'died_dt',
#'Died_Dt_str',
'earliest_cnty_dt',
#'Earliest_Cnty_Dt_str',
'earliest_state_dt',
#'Earliest_State_Dt_str',
#'Elect_Notif_CDC_Dt',
#'Elect_Notif_CDC_Dt_str',
#'Ethnicity',
'ethnicity_txt',
'first_elec_submit_dt',
#'First_Elec_Submit_Dt_str',
#'Generic_Version',
'hosp_admit_dt',
#'Hosp_Admit_Dt_str',
#'Hosp_Days',
#'Hosp_Dis_Dt',
#'Hosp_Dis_Dt_str',
#'Hospitalized',
#'Illness_Dur',
#'Illness_Dur_Units',
#'Illness_End_Dt',
#'Illness_End_Dt_str',
'illness_onset_dt',
#'Illness_Onset_Dt_str',
#'Immediate_NNC',
#'Import_City',
#'Import_City_txt',
#'Import_Cntry',
#'Import_Cntry_oth',
#'Import_Cnty',
#'Import_Cnty_txt',
#'Import_Code',
#'Import_Code_txt',
#'Import_State',
#'Import_State_txt',
'invest_start_dt',
#'Invest_Start_Dt_str',
#'Jurisdiction_Code',
#'Legacy_Case_ID',
#'Loc_Rec_ID',
#'Local_Subj_ID',
#'MMWR_Week',
#'MMWR_Year',
#'MMWR_Year_str',
#'Message_Con_ID',
#'Message_Dt',
#'Message_Dt_str',
#'Name_Report_CDC',
#'Name_Report_Email',
#'Nat_Report_Jurisdiction',
'notif_result_status',
#'Oth_Birth_Place',
#'Oth_Race',
#'Outbreak_Ind',
#'Outbreak_Name',
'phd_notif_dt',
#'PHD_Notif_Dt_Str',
#'Phone_Report_CDC',
'pregnant',
#'Proc_ID',
#'Profile_Version',
'report_county',
#'Report_County_txt',
'report_dt',
#'Report_Dt_str',
#'Report_Source_Type',
#'Report_Source_Type_oth',
#'Report_Source_Zip',
'report_state',
#'Report_State_txt',
#'Report_source_type_txt',
#'Sending_App_ID',
#'Sending_App_Name',
#'Sending_Fac_ID',
#'Sending_Fac_Name',
'sex',
#'State_Case_ID',
#'Transmission_Mode',
#'Transmission_Mode_oth',
#'Verbal_Notif_CDC_Dt',
#'Verbal_Notif_CDC_Dt_str',
#'birth_country_txt',
'current_record_flag',
#'import_cntry_txt',
#'message_seq_id',
#'msg_received_dttm',
'msg_transaction_id',
#'mvps_datetime_created',
#'mvps_datetime_updated',
#'nat_report_jurisdiction_txt',
#'transmission_mode_txt',
# TBD: 'Subj_County'
]


# In[5]:


name_map = {'age': 'Age',
            'age_units': 'Age_Units',
            'birth_date_str': 'Birth_Date_str',
            'case_status_txt': 'Case_Status_txt',
            'condition_code': 'Condition_Code',
            'diag_dt': 'Diag_Dt',
            'died_dt': 'Died_Dt',
            'earliest_cnty_dt': 'Earliest_Cnty_Dt',
            'earliest_state_dt': 'Earliest_State_Dt',
            'ethnicity_txt': 'Ethnicity_txt',
            'first_elec_submit_dt': 'First_Elec_Submit_Dt',
            'hosp_admit_dt': 'Hosp_Admit_Dt',
            'illness_onset_dt': 'Illness_Onset_Dt',
            'invest_start_dt': 'Invest_Start_Dt',
            'notif_result_status': 'Notif_Result_Status',
            'phd_notif_dt': 'PHD_Notif_Dt',
            'pregnant': 'Pregnant',
            'report_county': 'Report_County',
            'report_dt': 'Report_Dt',
            'report_state': 'Report_State',
            'sex': 'Sex'
}


# In[6]:


# File path constants
JURISDICTION_TO_STATE = '/data/csels/jurisdiction_to_state_mapping.csv'
OFFICIAL_FIPS = '/data/csels/us_official_fips_codes.csv'
NETSS_OUTPUT_ROOT = '/data/csels/preprocessed_data/netss_update'
HL7_OUTPUT_ROOT = '/data/csels/preprocessed_data/hl7_update'
NETSS_INPUT = '/data/csels/netss_oct_2020'
HL7_INPUT = '/data/csels/hl7_feb_2021'
GENERAL_INPUT = '/data/csels'


# In[7]:


# Define preferred column output orders and variables
netss_output_order = ['EVENTD','COUNT','AGE','AGETYPE','SEX','RACE','HISPANIC','CASSTAT','COUNTY']
hl7_output_order = [
    'Report_Dt_new','COUNT','Age','Age_Units','Sex','Ethnicity_txt',
    'race_mapped','Case_Status_txt', 'Birth_Date_str', 
    'Notif_Result_Status', 'Pregnant', 'Report_County', 
    'First_Elec_Submit_Dt','Subj_County', 'Diag_Dt', 'Died_Dt', 
    'Earliest_Cnty_Dt', 'Earliest_State_Dt', 'Hosp_Admit_Dt',
    'Illness_Onset_Dt', 'Invest_Start_Dt', 'PHD_Notif_Dt'
]


# In[8]:


# Load jurisdiction mappings
mapping = pd.read_csv(JURISDICTION_TO_STATE)
jurisdiction_code_to_name = dict(zip(mapping['Jurisdiction'], mapping['State']))


# In[9]:


#set of offical fips counties
fips = pd.read_csv(OFFICIAL_FIPS, dtype = str)
valid_county_set = set(fips['Concept Code'].values)


# ## Remove all existing preprocessed data
# This will remove all existing preprocessed data to give a clean slate for this run. It is done to protect against potential contamination that would result should this current run be aborted prematurely.

# In[10]:


# Remove the folders if they exist
if os.path.exists(NETSS_OUTPUT_ROOT):
    shutil.rmtree(NETSS_OUTPUT_ROOT)
    
if os.path.exists(HL7_OUTPUT_ROOT):
    shutil.rmtree(HL7_OUTPUT_ROOT)

# Make new folders
os.mkdir(NETSS_OUTPUT_ROOT)
os.mkdir(HL7_OUTPUT_ROOT)


# # Define preprocessing programs
# NETSS and HL7 data will be processed according to these function definitions.

# In[11]:


# function to apply to netss data
# modifies county codes based on conditions
def validate_county_code_netss(row):
    
    #identify row elements
    jur_num = row['STATE']
    county = row['COUNTY']
    
    #turn string to int
    if isinstance(county, str):
        county_num = int(county)
    #turn float to int except nan
    elif isinstance(county, float):
        if np.isnan(county):
            county_num = county
        else:
            county_num = int(county)
    else:
        county_num = county
    #return county if valid and between 0 and 999
    #else return 999
    if 0 < county_num < 999:
        county_try = f'{jur_num:02d}{county_num:03d}'
        if county_try in valid_county_set:
            return county_num
        else:
            return 999
    else:
        return 999


# In[12]:


def preprocess_netss_data(df):
    
    # Discard records where counts are not explicity set to 1, per CDC recommendations
    df = df.loc[df['COUNT'] == 1, :]
    
    # Discard records with missing EVENTD variables
    df = df.loc[np.logical_not(df['EVENTD'].isna()), :]
        
    # Update from 6/4/2020 - it was noted that EVENTD range from years 1899 to 2201.
    # We'll keep only complete years (2014-2019 for this currently being written in 2020). 
    # This is generalized such that if this script is run in 2021, 2020 will be included, and so on with subsequent years.
    valid_year = [x.year in range(2014,int(datetime.date.today().year)) for x in pd.to_datetime(df['EVENTD'], errors = 'coerce')]
    df = df.loc[valid_year, :]
    
    df['COUNTY'] = df.apply(validate_county_code_netss, axis = 1)
    
    # Obtain jurisdiction combinations
    jurisdictions = df['STATE'].unique()
    
    # Process each jurisdiction
    for j in jurisdictions:
        
        # Get the name of the jurisdiction
        jurisdiction_name = jurisdiction_code_to_name[int(j)]
        
        # Create a folder for this jurisdiction
        direc = (NETSS_OUTPUT_ROOT + '/{}').format(jurisdiction_name)
        os.mkdir(direc)
        
        # Subset to current jurisdiction
        df_jurisdiction = df.loc[df['STATE'] == j, :]
        
        # The STATE variable is no longer needed after subsetting
        df_jurisdiction = df_jurisdiction.drop(columns=['STATE'])
        
        # Obtain disease combinations for this jurisdiction
        diseases = df_jurisdiction['EVENT'].unique()
        
        # Process each disease within jurisdiction
        for d in diseases:
            
            # Subset to current disease
            df_disease = df_jurisdiction.loc[df['EVENT'] == d, :]
            
            # The EVENT variable is no longer needed after subsetting
            df_disease = df_disease.drop(columns=['EVENT'])
            
            # Initiate time series dataframe
            counts = df_disease.reset_index(drop = True)
    
            # Capture the time range and make daily increments for the range of data
            min_dt = pd.to_datetime(counts['EVENTD']).min()
            max_dt = pd.to_datetime(counts['EVENTD']).max()
            
            # Pad on the left and right to extend the series cleanly to the start/end of months in the data range
            series_min_dt = pd.Timestamp('{}-{}-01'.format(min_dt.year,min_dt.month))
            series_max_dt = max_dt + pd.tseries.offsets.MonthEnd(0)
            all_dates = pd.DataFrame(pd.date_range(start = series_min_dt, end = series_max_dt), columns = ['EVENTD'])
            
            # Merge the counts to the padded data and set missing counts to 0
            counts['EVENTD'] = pd.to_datetime(counts['EVENTD'])
            time_series = all_dates.merge(counts, on = ['EVENTD'], how = 'left')
            time_series['COUNT'] = time_series['COUNT'].fillna(0)
            
            # Export this dataset as a csv file within the jurisdiction directory
            filename = '{}.csv'.format(d)
            filepath = os.path.join(NETSS_OUTPUT_ROOT, jurisdiction_name, filename)
            time_series[netss_output_order].to_csv(filepath, index = False)
            #time_series[netss_output_order].to_csv('./netss_update/{}/{}.csv'.format(jurisdiction_name, d), index = False)


# In[13]:


#function to apply to hl7 data
#modifies county codes based on condition
def validate_county_code_hl7_rep(row):
    
    #identify row elements
    jur_num = row['Report_State']
    county = row['Report_County']
    
    #try if data is string else return '99999'
    try:
        #must be five digit string
        if re.fullmatch(r'\d{5}', county):
            #county jurisdiction
            jur_test = county[:2]
            #check if county jurisdiction matches actual
            if jur_test == f'{jur_num:02d}':
                #keep county if valid
                if county in valid_county_set:
                    return county
                #keep state jurisdiction if match
                else:
                    return f'{jur_num:02d}999'
            #jurisdiction doesn't match
            else:
                return jur_test + '999'
        #replace non five digit string
        else:
            return '99999'
    except:
        return '99999'


# In[14]:


def validate_county_code_hl7_subj(row):
    
    #identify row elements
    jur_num = row['Report_State']
    county = row['Subj_County']
    
    #try if data is string else return '99999'
    try:
        #must be five digit string
        if re.fullmatch(r'\d{5}', county):
            #county jurisdiction
            jur_test = county[:2]
            #check if county jurisdiction matches actual
            if jur_test == f'{jur_num:02d}':
                #keep county if valid
                if county in valid_county_set:
                    return county
                #keep state jurisdiction if match
                else:
                    return f'{jur_num:02d}999'
            #jurisdiction doesn't match
            else:
                return jur_test + '999'
        #replace non five digit string
        else:
            return '99999'
    except:
        return '99999'


# In[15]:


def preprocess_hl7_data(df):
    
    # Exclude rows where current_record_flag is not 'Y'
    df = df.loc[df['current_record_flag'] == 'Y', :]
    df = df.drop(columns = 'current_record_flag')
    
    # Exclude rows where report_dt is missing
    df = df.copy() # to avoid SettingWithCopyWarning
    lookup = df.loc[:, ['Report_Dt', 'PHD_Notif_Dt', 'Earliest_State_Dt', 'Earliest_Cnty_Dt']].notnull().idxmax(1)
    idx, cols = pd.factorize(lookup)
    df['Report_Dt_new'] = df.reindex(cols, axis = 1).to_numpy()[np.arange(len(df)), idx]
    #df = df.loc[np.logical_not(df['Report_Dt_new'].isna()), :]
    
    # Calculate the total number of observations for each year and exclude if < 1000
    report_dt = pd.to_datetime(df['First_Elec_Submit_Dt'])
    report_yr = [x.year for x in report_dt]
    cnt = collections.Counter(report_yr)
    valid_years = [x for x, count in cnt.items() if count >= 1000]
    df = df.loc[pd.Series(report_yr).isin(valid_years).values, :]
    
    # Replace report_dt with the date part of the full string expression
    df['First_Elec_Submit_Dt'] = df['First_Elec_Submit_Dt'].str.slice(stop=10)
    
    # Set COUNT variable to 1 (to align with the NETSS 1 row = 1 COUNT methodology)
    df['COUNT'] = 1
    
    df = df[df['Report_State'].notna()]
    df['Report_State'] = df['Report_State'].astype(int)
    
    df['Report_County'] = df.apply(validate_county_code_hl7_rep, axis = 1)
    df['Subj_County'] = df.apply(validate_county_code_hl7_subj, axis = 1)
    
    # Obtain jurisdiction combinations
    jurisdictions = df['Report_State'].unique()
    
    # Process each jurisdiction
    for j in jurisdictions:
        
        # Get the name of the jurisdiction
        jurisdiction_name = jurisdiction_code_to_name[int(j)]
        
        # Create a folder for this jurisdiction
        direc = (HL7_OUTPUT_ROOT + '/{}').format(jurisdiction_name)
        os.mkdir(direc)
        
        # Subset to current jurisdiction
        df_jurisdiction = df.loc[df['Report_State'] == j, :]
        
        # The report_state variable is no longer needed after subsetting
        df_jurisdiction = df_jurisdiction.drop(columns=['Report_State'])
        
        # Obtain disease combinations for this jurisdiction
        diseases = df_jurisdiction['Condition_Code'].unique()
        
        # Process each disease within jurisdiction
        for d in diseases:
            
            # Subset to current disease
            df_disease = df_jurisdiction.loc[df['Condition_Code'] == d, :]
            
            # The condition_code variable is no longer needed after subsetting
            df_disease = df_disease.drop(columns=['Condition_Code'])
            
            # Initiate time series dataframe
            counts = df_disease.reset_index(drop = True)
    
            # Capture the time range and make daily increments for the range of data
            min_dt = pd.to_datetime(counts['First_Elec_Submit_Dt']).min()
            max_dt = pd.to_datetime(counts['First_Elec_Submit_Dt']).max()
            
            # Pad on the left and right to extend the series cleanly to the start/end of months in the data range
            series_min_dt = pd.Timestamp('{}-{}-01'.format(min_dt.year,min_dt.month))
            series_max_dt = max_dt + pd.tseries.offsets.MonthEnd(0)
            all_dates = pd.DataFrame(pd.date_range(start = series_min_dt, end = series_max_dt), 
                                     columns = ['First_Elec_Submit_Dt'])
            
            # Merge the counts to the padded data and set missing counts to 0
            counts['First_Elec_Submit_Dt'] = pd.to_datetime(counts['First_Elec_Submit_Dt'])
            time_series = all_dates.merge(counts, on = ['First_Elec_Submit_Dt'], how = 'left')
            time_series['COUNT'] = time_series['COUNT'].fillna(0)
            
            # Export this dataset as a csv file within the jurisdiction directory
            filename = '{}.csv'.format(d)
            filepath = os.path.join(HL7_OUTPUT_ROOT, jurisdiction_name, filename)
            time_series[hl7_output_order].to_csv(filepath, index = False)
            #time_series[hl7_output_order].to_csv('./hl7_update/{}/{}.csv'.format(jurisdiction_name, d), index = False)    


# In[16]:


# Date parsing utility function for the function to follow
def dateParse(df, x):
    
    # It seems that illness_end_dt_str are formatted as YYYYMMDD, while other datetime variables are YYYY-MM-DD.
    # However, there are other peculiaries, such as illness_end_dt_str ending in .000 and datetime variables
    # requiring their time component to be stripped away, which is hanlded by the string slicing below.
    if x == 'Illness_End_Dt_str':
        return pd.to_datetime(df[x].str.slice(start = 0, stop = 8), errors='coerce', format = '%Y%m%d')
    else:
        return pd.to_datetime(df[x].str.slice(start = 0, stop = 10), errors='coerce', format = '%Y-%m-%d')

# Function to take nonmissing dates in reference to their days from Report_Dt, store as a tuple, and count unique values    
def create_tuple_csv(df):
    
    # Exclude rows where report_dt is missing
    df = df.copy() # to avoid SettingWithCopyWarning
    lookup = df.loc[:, ['Report_Dt', 'PHD_Notif_Dt', 'Earliest_State_Dt', 'Earliest_Cnty_Dt']].notnull().idxmax(1)
    idx, cols = pd.factorize(lookup)
    df['Report_Dt_new'] = df.reindex(cols, axis = 1).to_numpy()[np.arange(len(df)), idx]
    #df = df.loc[np.logical_not(df['Report_Dt_new'].isna()), :]
    
    # Necessary variables
    core_dts = df[['Diag_Dt', 'Died_Dt', 'First_Elec_Submit_Dt', 'Hosp_Admit_Dt', 'Illness_Onset_Dt', 'Invest_Start_Dt',
                  'Report_Dt_new', 'current_record_flag']]
    
    # Apply restrictions that we intend to apply before generating synthetic data

    # Exclude rows where current_record_flag is not 'Y'
    core_dts = core_dts.loc[core_dts['current_record_flag'] == 'Y', :]
    core_dts = core_dts.drop(columns = 'current_record_flag')

    # Exclude rows where report_dt is missing
    core_dts = core_dts.copy() # to avoid SettingWithCopyWarning
    #core_dts = core_dts.loc[np.logical_not(core_dts['Report_Dt_new'].isna()), :]
    
    # Parse the dates to string format
    for i in ['Diag_Dt', 'Died_Dt', 'First_Elec_Submit_Dt', 'Hosp_Admit_Dt', 'Illness_Onset_Dt', 'Invest_Start_Dt', 
              'Report_Dt_new']:
        core_dts[i] = dateParse(core_dts, i)
    
    # List of each individual date difference tuple 
    # Only tuples where all values are nonnegative are kept
    tupes = []
    for index, row in core_dts.iterrows():
        report = row['First_Elec_Submit_Dt']
        if row['Diag_Dt'] != 'None': 
            a = (row['Diag_Dt'] - report) / np.timedelta64(1, 'D')
        else: 
            a = 'None'
        if row['Died_Dt'] != 'None': 
            b = (row['Died_Dt'] - report) / np.timedelta64(1, 'D') 
        else: 
            b = 'None'
        if row['Report_Dt_new'] != 'None': 
            c = (row['Report_Dt_new'] - report) / np.timedelta64(1, 'D') 
        else: 
            c = 'None'
        if row['Hosp_Admit_Dt'] != 'None': 
            d = (row['Hosp_Admit_Dt'] - report) / np.timedelta64(1, 'D') 
        else: 
            d = 'None'
        if row['Illness_Onset_Dt'] != 'None': 
            e = (row['Illness_Onset_Dt'] - report) / np.timedelta64(1, 'D') 
        else: 
            e = 'None'
        if row['Invest_Start_Dt'] != 'None': 
            f = (row['Invest_Start_Dt'] - report) / np.timedelta64(1, 'D') 
        else: 
            f = 'None'
        tupes.append((a,b,c,d,e,f))
    
    # Counter object for each individual tuple           
    counts = collections.Counter(tupes)
    hist = pd.DataFrame.from_dict(counts, orient='index').reset_index()
    hist[['d1','d2','d3','d4','d5','d6']] = pd.DataFrame(hist['index'].tolist(), index=hist.index)
    hist = hist.drop('index', axis = 1)
    hist = hist.rename(columns = {0: 'count'})
    
    # Csv file writer
    filename = 'tuple_histogram.csv'
    filepath = os.path.join(HL7_OUTPUT_ROOT, filename)
    hist.to_csv(filepath)
    #hist.to_csv('./hl7_update/tuple_histogram.csv')


# We also define a utility function to transform the `gen_race.csv` file into a condensed (1 row = 1 `msg_transaction_id`) format. This allows it to be joined to core, which adds a single column called `race_mapped`. This column is a concatenation of all races given by the corresponding `msg_transaction_id` in the race table.
# 
# Please refer to this reference when considering the race code mappings:
# https://phinvads.cdc.gov/vads/ViewValueSet.action?oid=2.16.840.1.114222.4.11.7205

# In[17]:


def make_condensed_race_table():
    
    # Load race table
    filename = 'genRace.csv'
    filepath = os.path.join(HL7_INPUT, filename)
    race = pd.read_csv(filepath, dtype = str, encoding = 'cp1252', quotechar = '^')
    #race = pd.read_csv('../hl7_feb_2021/genRace.csv',
    #                   dtype = str, encoding = 'cp1252', quotechar = '^')

    
    # Define race code mapping to get codes aligned with standard concepts
    race_code_map = {
        '2106-3':'2106-3',               # White
        'UNK':'NullFlavor',              # NullFlavor (Unknown)
        '2054-5':'2054-5',               # Black or African American
        '2131-1':'2131-1',               # Other Race
        '2028-9':'2028-9',               # Asian
        '1002-5':'1002-5',               # American Indian or Alaska Native
        '2076-8':'2076-8',               # Native Hawaiian or Other Pacific Islander
        'C':'NullFlavor',                # NullFlavor (Unknown)
        'PHC1175':'NullFlavor',          # NullFlavor (Refused to answer)
        'U':'NullFlavor',                # NullFlavor (Unknown)
        'NASK':'NullFlavor',             # NullFlavor (Not Asked)
        'H':'NullFlavor',                # NullFlavor (Unknown)
        'M':'NullFlavor',                # NullFlavor (Unknown)
        'White':'2106-3',                # White
        'X':'NullFlavor',                # NullFlavor (Unknown)
        'R':'NullFlavor',                # NullFlavor (Unknown)
        '2054-4':'2054-5',               # Black or African American (Presumed inteded to be 2054-5)
        'WHITE (CAUCASIAN)':'2106-3',    # White
        'WH':'2106-3',                   # White
        'U':'NullFlavor',                # NullFlavor (Presumed Unknown)
        'AI':'1002-5',                   # Presumed AI stands for American Indian
        'NH':'2076-8',                   # Presumed NH stands for Native Hawaiian
        'NOT PROVIDED':'NullFlavor',     # NullFlavor
        'UNKNOWN':'NullFlavor',          # NullFlavor
        'Unknown':'NullFlavor',          # NullFlavor
        'PHC1367':'NullFlavor',          # NullFlavor (Refused)
        'P':'NullFlavor',                # NullFlavor (Unknown)
        'Pt Declined':'NullFlavor',      # NullFlavor
        'CAUCASIAN':'2106-3',            # White
        'CA':'2106-3',                   # White
        'Black':'2054-5',                # Black or African American
        'N':'NullFlavor',                # NullFlavor
        'Not Reported':'NullFlavor',     # NullFlavor
        '2058-6':'2131-1',               # Other (Undocumented code)
        'OTHER NONWHITE':'2131-1',       # Other
        'Black or African Ame':'2054-5', # Black or African American
        '1004-1':'1002-5',               # American Indian
        'ASKU':'NullFlavor',             # NullFlavor (Asked but unknown)
        '1':'NullFlavor',                # NullFlavor (Unknown)
        '2186-5':'2131-1',               # Other (Undocumented code)
        '2076-8s':'2131-1',              # Other (Undocumented code)
        'code failure':'NullFlavor',     # NullFlavor
        'AFRICAN AMERICAN OR':'2054-5',  # Black or African American
        'Other Race':'2131-1',           # Other
        'Hispanic':'2131-1',             # Other
        'Hispanic or Latino':'2131-1',   # Other
        'FALSE':'NullFlavor',            # NullFlavor
        'ASIAN':'2028-9',                 # Asian
        #new
        '1002-5, 2054-5, 2106': '2131-1',
        '2029-7': '2028-9',
        '2036-2': '2028-9',
        '2039-6': '2028-9',
        '2086-7': '2028-9',
        'AA': '2054-5',
        'AS': '2028-9',
        'Asian': '2028-9',
        'CAC': 'NullFlavor',
        'CACH': 'NullFlavor',
        'CAF': 'NullFlavor',
        'CAG': 'NullFlavor',
        'CAH': 'NullFlavor',  
        'CAI': 'NullFlavor',
        'CAJ': 'NullFlavor',
        'CAK': 'NullFlavor',
        'CAL': 'NullFlavor',
        'CAM': 'NullFlavor',
        'CAMI': 'NullFlavor',
        'CAMP': 'NullFlavor',
        'CAN': 'NullFlavor',
        'CAO': 'NullFlavor',
        'CAOT': 'NullFlavor',
        'CAP': 'NullFlavor',
        'CAS': 'NullFlavor',
        'CAV': 'NullFlavor',
        'Caucasian': '2106-3',
        'FALSE': 'NullFlavor',
        'FEMALE': 'NullFlavor',
        'O': '2131-1',
        'OTHER ASIAN': '2131-1',
        'Other': '2131-1',
        'PHC1367': 'NullFlavor',
        'S': 'NullFlavor',
        'WHITE': '2016-3',
        ' ': 'NullFlavor'
    }

    # Apply mapping
    race['race_mapped'] = race['race'].map(race_code_map)

    # Calculate condensed race by group
    g = race.groupby('msg_transaction_id')['race_mapped']
    race_condensed = g.apply(lambda x: ';'.join(sorted(x.unique()))).reset_index()
    
    # A subsequent step is required as an artifact of the NullFlavor mappings
    # Since an individual can be marked, for example, as both "unknown" and "refused", they can have two or more NullFlavor entries when condensed
    # We will map such cases to a single NullFlavor entry.
    is_nullflavor = race_condensed['race_mapped'].map(lambda x: set(x.split(';')) == {'NullFlavor'})
    race_condensed.loc[is_nullflavor, 'race_mapped'] = 'NullFlavor'
    
    # We will also condense string sequences of the form x1;x2;...;xn;NullFlavor to be x1;x2;...;xn
    # If at least one race is indicated, then we posit that a NullFlavor suffix is not necessary
    race_condensed['race_mapped'] = race_condensed['race_mapped'].map(lambda x: ';'.join([item for item in x.split(';') if item not in 'NullFlavor']) if ((len(x.split(';')) > 1) and ('NullFlavor' in x.split(';'))) else x)

    # This dataset is now ready to be joined to core
    return race_condensed


# ## Load, preview, and process data

# ### NETSS

# In[19]:


file1 = os.path.join(NETSS_INPUT, 'year2015a.csv')
file2 = os.path.join(NETSS_INPUT, 'year2016a.csv')
file3 = os.path.join(NETSS_INPUT, 'year2017pt1.csv')
file4 = os.path.join(NETSS_INPUT, 'year2017pt2a.csv')
file5 = os.path.join(NETSS_INPUT, 'year2018.csv')
file6 = os.path.join(NETSS_INPUT, 'nndss19.csv')
file7 = os.path.join(GENERAL_INPUT, 'nndss20_update.csv')
netss1 = pd.read_csv(file1, usecols = usecols_NETSS,
                    dtype = {'EVENTD':str,'COUNT':int,'STATE':int,'EVENT':int})
netss2 = pd.read_csv(file2, usecols = usecols_NETSS,
                    dtype = {'EVENTD':str,'COUNT':int,'STATE':int,'EVENT':int})
netss3 = pd.read_csv(file3, usecols = usecols_NETSS,
                    dtype = {'EVENTD':str,'COUNT':int,'STATE':int,'EVENT':int})
netss4 = pd.read_csv(file4, usecols = usecols_NETSS,
                    dtype = {'EVENTD':str,'COUNT':int,'STATE':int,'EVENT':int})
netss5 = pd.read_csv(file5, usecols = usecols_NETSS,
                    dtype = {'EVENTD':str,'COUNT':int,'STATE':int,'EVENT':int})
netss6= pd.read_csv(file6, usecols = usecols_NETSS,
                    dtype = {'EVENTD':str,'COUNT':int,'STATE':int,'EVENT':int})
netss7 = pd.read_csv(file7, usecols = usecols_NETSS,
                    dtype = {'EVENTD':str,'COUNT':int,'STATE':int,'EVENT':int})
netss = pd.concat([netss1,netss2,netss3,netss4,netss5,netss6,netss7], ignore_index = True)


# In[20]:


# Preview NETSS data head (Optional)
# netss.head(10)


# In[21]:


# Preview NETSS data tail (Optional)
# netss.tail(10)


# In[22]:


# Process netss data
preprocess_netss_data(netss)


# ### HL7 - gen_core

# In[18]:


# Load data subset to relevant columns for processing
#hl7 = pd.read_csv('../hl7_feb_2021/genCore_Export_y2.zip', usecols = usecols_HL7,
                  #quotechar='^', dtype = str, header = (0), encoding = 'cp1252', compression = 'zip')
file1 = os.path.join(HL7_INPUT, 'genCore_2019_Export_transid.csv')
file2 = os.path.join(HL7_INPUT, 'genCore_2020_Export_transid.csv')
hl719 = pd.read_csv(file1, header = (1), quotechar = '^', dtype = str, 
                    encoding = 'cp1252', usecols = usecols_HL7)
hl720 = pd.read_csv(file2, quotechar = '^', dtype = str, encoding = 'cp1252',
                    usecols = usecols_HL7)
hl7 = pd.concat([hl719,hl720], ignore_index = True)
hl7 = hl7.rename(mapper = name_map, axis = 1)
#hl7 = pd.read_csv('../hl7_oct_2020/GenCore.csv',
#                    usecols = usecols_HL7,
#                    dtype = str)


# In[19]:


# Create the condensed race table -- If this fails with a TypeError, it implies that not all of the race categories were accounted for
try:
    race_condensed = make_condensed_race_table()
except TypeError:
    print("Please ensure that 'race_code_map' contains all of the mappings below:")
    file = os.path.join(HL7_INPUT, 'genRace.csv')
    df_race = pd.read_csv(file, usecols = ['race'], dtype = str, encoding = 'cp1252', quotechar = '^')
    for x in sorted(df_race['race'].unique()):
        print(x)


# In[20]:


# Join this to core
hl7 = hl7.merge(race_condensed, on = 'msg_transaction_id', how = 'left')


# In[21]:


# Read in file to merge with hl7
# Only need one variable
file = os.path.join(HL7_INPUT, 'genAddress.csv')
subj_county = pd.read_csv(file, 
usecols = ['msg_transaction_id', 'Subj_County'], 
dtype = str, encoding = 'cp1252', quotechar = '^').dropna().reset_index(drop = True)

# Joining the variable to the core
hl7 = hl7.merge(subj_county, on = 'msg_transaction_id', how = 'left')


# In[22]:


# A subsequent step here is required to address cases where the msg_transaction_id in core was absent from the race table.
# Such cases will be mapped to NullFlavor.
hl7.loc[hl7['race_mapped'].isna(),'race_mapped'] = 'NullFlavor'


# In[23]:


# Preview HL7 data head (Optional)
# hl7.head(10)


# In[24]:


# Preview HL7 data tail (Optional)
# hl7.tail(10)


# In[25]:


# Process HL7 data
preprocess_hl7_data(hl7)


# In[26]:


# Create the tuple histogram csv
#create_tuple_csv(hl7)


# In[27]:


# End timer
elapsed_time = time.time() - start_time
print("Time elapsed: ", round(elapsed_time/3600,2), "Hours")


# In[ ]:




