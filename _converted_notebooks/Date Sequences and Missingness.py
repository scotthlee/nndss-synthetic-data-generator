#!/usr/bin/env python
# coding: utf-8

# ## Date Sequences and Missingness
# Per a CDC/GTRI meeting on 10/22/2020, we have decided to examine the HL7 data dates more closely for the purpose of synthetic data generation. Generation around dates require special considerations because certain logical orders/relationships may need to be preserved when generating them synthetically (e.g. One can't be discharged from the hospital before being admitted). 
# 
# The notebook explores the following:
# - We need to understand missingness overall and by group (by jurisdiction and by condition).
# - We also need to understand whether dates can appear before/after one another in the data.

# In[1]:


# Imports
import collections
import numpy as np
import pandas as pd
import seaborn as sns
from itertools import combinations


# In[2]:


# Options to display full tables instead of ellipses
pd.set_option("display.max_columns",50)
pd.set_option("display.max_rows",5000)


# In[3]:


# Set path to HL7 Gen Core dataset
GEN_CORE = '../hl7_oct_2020/GenCore.csv'


# In[4]:


# These are all of the variables which appear to be a date in the Hl7 data
date_vars = ['Birth_Date',
             'Illness_Onset_Dt',
             'Illness_End_Dt_str',
             'Diag_Dt',
             'Hosp_Admit_Dt',
             'Hosp_Dis_Dt',
             'Invest_Start_Dt',
             'First_Elec_Submit_Dt',
             'Elect_Notif_CDC_Dt',
             'Report_Dt',
             'Earliest_Cnty_Dt',
             'Earliest_State_Dt',
             'MMWR_Year',
             'PHD_Notif_Dt',
             'Message_Dt',
             'msg_received_dttm',
             'mvps_datetime_created',
             'mvps_datetime_updated',
             'current_record_flag', # Here and below are not dates but still needed for groupby operations
             'Report_State_txt',
             'Condition_Code_txt']


# In[5]:


core_dts = pd.read_csv(GEN_CORE, usecols = date_vars, dtype = str)


# In[6]:


# Apply restrictions that we intend to apply before generating synthetic data

# Exclude rows where current_record_flag is not 'Y'
core_dts = core_dts.loc[core_dts['current_record_flag'] == 'Y  ', :]
core_dts = core_dts.drop(columns = 'current_record_flag')

# Exclude rows where report_dt is missing
#core_dts = core_dts.copy() # to avoid SettingWithCopyWarning
core_dts = core_dts.loc[np.logical_not(core_dts['Report_Dt'].isna()), :]


# ### Overall Missingness

# In[7]:


# Missingness "Heatmap" -- Missing is ivory-colored, completed is black
sns.set(rc={'figure.figsize':(16,32)})
sns.heatmap(core_dts.isnull(), cbar=False).figure.savefig("./figures/overall_missingness.png")


# In[8]:


# Missing percentages by variable
core_dts.apply(lambda x: sum(x.isnull()) / core_dts.shape[0] * 100)


# In[9]:


# So, we can set aside these dates that are never missing (or very rarely missing)
not_missing = ['First_Elec_Submit_Dt',
               'Elect_Notif_CDC_Dt',
               'Report_Dt',
               'MMWR_Year',
               'Message_Dt',
               'msg_received_dttm',
               'mvps_datetime_created']


# In[10]:


missing = list(set(date_vars) - set(not_missing) - set(['current_record_flag', 'Report_State_txt', 'Condition_Code_txt']))
missing


# In[11]:


# Define groupbys

# By diagnosis
dx = core_dts.groupby('Condition_Code_txt')

# By jurisdiction
dj = core_dts.groupby('Report_State_txt')

# By both
dx_dj = core_dts.groupby(['Report_State_txt', 'Condition_Code_txt'])
dj_dx = core_dts.groupby(['Condition_Code_txt', 'Report_State_txt'])


# ### Missingness by Condition

# In[12]:


# Missingness by condition
missing_by_condition = pd.concat([dx[y].apply(lambda x: sum(x.isnull()) / len(x) * 100) for y in missing], axis = 1)
missing_by_condition['sample_size'] = dx.size().values
missing_by_condition


# In[13]:


# Export table
missing_by_condition.to_csv('./datasets/missing_by_condition.csv')


# In[14]:


sns.set(rc={'figure.figsize':(16,48)})
sns.heatmap(missing_by_condition.drop(columns='sample_size'), cbar = False, cmap = "mako").figure.savefig("./figures/missing_by_condition.png")


# ### Missingness by Jurisdiction

# In[15]:


# Missingness by jurisdiction
missing_by_jurisdiction = pd.concat([dj[y].apply(lambda x: sum(x.isnull()) / len(x) * 100) for y in missing], axis = 1)
missing_by_jurisdiction['sample_size'] = dj.size().values
missing_by_jurisdiction


# In[16]:


# Export table
missing_by_jurisdiction.to_csv('./datasets/missing_by_jurisdiction.csv')


# In[17]:


sns.set(rc={'figure.figsize':(16,20)})
sns.heatmap(missing_by_jurisdiction.drop(columns='sample_size'), cbar = False, cmap = "mako").figure.savefig("./figures/missing_by_jurisdiction.png")


# ### Missingness by Condition and Jurisdiction

# In[18]:


missing_by_condition_and_jurisdiction = pd.concat([dx_dj[y].apply(lambda x: sum(x.isnull()) / len(x) * 100) for y in missing], axis = 1)
missing_by_condition_and_jurisdiction['sample_size'] = dx_dj.size().values
missing_by_condition_and_jurisdiction


# In[19]:


# Export table
missing_by_condition_and_jurisdiction.to_csv('./datasets/missing_by_jurisdiction.csv')


# In[20]:


sns.set(rc={'figure.figsize':(16,200)})
sns.heatmap(missing_by_condition_and_jurisdiction.drop(columns='sample_size'), cbar = False, cmap = "mako").figure.savefig("./figures/missing_by_condition_jurisdiction.png")


# ### Missingness by Jurisdiction and Condition

# In[21]:


missing_by_jurisdiction_and_condition = pd.concat([dj_dx[y].apply(lambda x: sum(x.isnull()) / len(x) * 100) for y in missing], axis = 1)
missing_by_jurisdiction_and_condition['sample_size'] = dj_dx.size().values
missing_by_jurisdiction_and_condition


# In[22]:


# Export table
missing_by_jurisdiction_and_condition.to_csv('./datasets/missing_by_jurisdiction_and_condition.csv')


# In[23]:


sns.set(rc={'figure.figsize':(16,200)})
sns.heatmap(missing_by_jurisdiction_and_condition.drop(columns='sample_size'), cbar = False, cmap = "mako").figure.savefig("./figures/missing_by_jurisdiction_condition.png")


# ## Date Sequences
# For any two date variable tuples (a, b), whenever both a and b are present, is a < b?

# In[24]:


# Gather date tuples
date_tuples = list(combinations(list(set(date_vars) - set(['current_record_flag', 'Report_State_txt', 'Condition_Code_txt'])), 2))
date_tuples


# In[25]:


# Date parsing utility function for the function to follow
def dateParse(df, x):
    
    # It seems that illness_end_dt_str are formatted as YYYYMMDD, while other datetime variables are YYYY-MM-DD.
    # However, there are other peculiaries, such as illness_end_dt_str ending in .000 and datetime variables
    # requiring their time component to be stripped away, which is hanlded by the string slicing below.
    if x == 'Illness_End_Dt_str':
        return pd.to_datetime(df[x].str.slice(start = 0, stop = 8), errors='coerce', format = '%Y%m%d')
    else:
        return pd.to_datetime(df[x].str.slice(start = 0, stop = 10), errors='coerce', format = '%Y-%m-%d')


# In[26]:


# Function to process a date tuple and assess how often both dates were present, and if they were, how often a < b
def processTuple(tup):
    
    # Unpack
    a, b = tup
       
    # Subset the data to complete cases (both a and b exist)
    sub = core_dts[[a,b]] 
    sub_complete = sub.loc[~(sub[a].isna() | sub[b].isna()), :]
    
    # Capture the sample size of this and the percentage of the data this represents
    complete_cases = sub_complete.shape[0]
    complete_case_pct = round(complete_cases / sub.shape[0] * 100,2)
    
    # Depending on the variable type, it needs to be converted to a date differently because the dates are in varying formats
    a_dates = dateParse(sub_complete, a)
    b_dates = dateParse(sub_complete, b)
    
    # Subtract and convert date differences to days
    a_b = (a_dates - b_dates) / np.timedelta64(1, 'D')
    
    # Gather statistics about this difference
    a_b_mean   = a_b.mean()
    a_b_std    = a_b.std()
    a_b_median = a_b.median()
    a_b_min    = a_b.min()
    a_b_max    = a_b.max()
    a_b_l95    = a_b.quantile(0.025)
    a_b_u95    = a_b.quantile(0.975)
    
    # Pass back a dataframe containing all difference information
    return pd.DataFrame(data = {'Date1':[a],
                                'Date2':[b],
                                'Complete_Cases':[complete_cases],
                                'Complete_Case_Pct':[complete_case_pct],
                                'Min':[a_b_min],
                                'Max':[a_b_max],
                                'Mean':[a_b_mean],
                                'Std':[a_b_std],
                                'Median':[a_b_median],
                                'L95':[a_b_l95],
                                'U95':[a_b_u95]})
    


# In[27]:


# Example of use of this function on a single tuple
processTuple(('Illness_Onset_Dt','Hosp_Dis_Dt'))


# In[28]:


# Call this on all date combination tuples
date_sequences = pd.concat([processTuple(x) for x in date_tuples])


# In[29]:


# Show the distributions of pairwise date differences (Date1 - Date2)
# i.e. NEGATIVE differences imply Date1 came BEFORE Date2
#      POSITIVE differences imply Date1 came AFTER  Date2
date_sequences


# In[30]:


# Export table
date_sequences.to_csv('./datasets/date_sequence_distributions.csv', index = False)

