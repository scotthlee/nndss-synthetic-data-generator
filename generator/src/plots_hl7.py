import os
import sys
import datetime
import numpy as np
import matplotlib.pyplot as plt

from src import hl7 as HL7

_SEX_TICK_LOCS = [0, 0.5, 1, 1.5, 2, 2.5, 3]
_SEX_TICK_LABELS = ('', 'U', '', 'M', '', 'F', '')
_SEX_ECDF_TICK_LABELS = ('U (0)', '', 'M (1)', '', 'F (2)', '', '')

_ETH_TICK_LOCS = [0,  0.5, 1,  1.5, 2,  2.5, 3,  3.5, 4]
_ETH_TICK_LABELS = ('', 'U', '', 'Y', '', 'N', '', 'O', '')
_ETH_ECDF_TICK_LABELS = ('U (0)', '', 'Y (1)', '', 'N (2)', '', 'O (3)', '', '')

_CS_TICK_LOCS = [0,  0.5, 1,  1.5, 2,  2.5, 3,  3.5, 4,  4.5, 5]
_CS_TICK_LABELS = ('', 'U', '', 'N', '', 'S', '', 'P', '', 'C', '')
_CS_ECDF_TICK_LABELS  = ('U (0)', '', 'N (1)', '', 'S (2)', '', 'P (3)', '', 'C (4)', '', '')

_PREG_TICK_LOCS = [0,  0.5, 1,  1.5, 2,  2.5, 3]
_PREG_TICK_LABELS = ('', 'U', '', 'N', '', 'Y', '')
_PREG_ECDF_TICK_LABELS = ('U (0)', '', 'N (1)', '', 'Y (2)', '', '')


###############################################################################
def plot_marginals(state):

    race_max = state.race_max
    county_max = state.county_max

    age_data = state.age_data
    sex_data = state.sex_data
    race_data = state.race_data
    ethnicity_data = state.ethnicity_data
    county_data = state.county_data
    case_status_data = state.case_status_data
    pregnant_data = state.pregnant_data

    AGE_MIN  = HL7.AGE_MIN
    AGE_MAX  = HL7.AGE_MAX
    SEX_MIN  = HL7.SEX_MIN
    SEX_MAX  = HL7.SEX_MAX
    RACE_MIN = HL7.RACE_MIN
    ETHNICITY_MIN = HL7.ETHNICITY_MIN
    ETHNICITY_MAX = HL7.ETHNICITY_MAX
    CASE_STATUS_MIN = HL7.CASE_STATUS_MIN
    CASE_STATUS_MAX = HL7.CASE_STATUS_MAX
    COUNTY_MIN = HL7.COUNTY_MIN
    PREGNANT_MIN = HL7.PREGNANT_MIN
    PREGNANT_MAX = HL7.PREGNANT_MAX

    # race_max is determined from the data, not from a file of constants
    RACE_MAX = race_max

    # county_max is determined from the data, not from a file of constants
    COUNTY_MAX = county_max

    # width, height    
    plt.figure(figsize=(20, 10))
    plt.hist(age_data, bins=AGE_MAX-AGE_MIN+1, range=(AGE_MIN, AGE_MAX+1),
             density=True, align='left', edgecolor='white')
    plt.xlabel('Age in Years', fontsize=12)
    plt.ylabel('Probability', fontsize=12)
    plt.title('AGE Marginal PDF', fontsize=16)
    plt.grid(zorder=-5)

    plt.figure(figsize=(20, 10))
    if COUNTY_MAX > 80:
        # don't use a white edgecolor; bins can be erased for large bin numbers
        plt.hist(county_data, bins=COUNTY_MAX-COUNTY_MIN+1,
                 range=(COUNTY_MIN, COUNTY_MAX+1), density=True)
    else:
        plt.hist(county_data, bins=COUNTY_MAX-COUNTY_MIN+1,
                 range=(COUNTY_MIN, COUNTY_MAX+1), density=True, edgecolor='white')
    plt.xlabel('Remapped County Code', fontsize=12)
    plt.ylabel('Probability', fontsize=12)
    plt.title('COUNTY Marginal PDF', fontsize=16)
    plt.grid(zorder=-5)

    plt.figure(figsize=(20, 10))
    if RACE_MAX > 80:
        plt.hist(race_data, bins=RACE_MAX-RACE_MIN+1,
                 range=(RACE_MIN, RACE_MAX+1), density=True)
    else:
        plt.hist(race_data, bins=RACE_MAX-RACE_MIN+1,
                 range=(RACE_MIN, RACE_MAX+1), density=True, edgecolor='white')
    if RACE_MAX < 20:
        # one tick for every value
        plt.xticks(range(0, RACE_MAX+1, 1))
    elif RACE_MAX < 40:
        # one tick every two values
        plt.xticks(range(0, RACE_MAX+1, 2))
    else:
        # use the default
        pass
    plt.xlabel('Remapped Race Code', fontsize=12)
    plt.ylabel('Probability', fontsize=12)
    plt.title('RACE Marginal PDF', fontsize=16)
    plt.grid(zorder=-5)

    # smaller plots
    plt.figure(figsize=(20,6))
    plt.subplot(1,4,1)
    plt.hist(sex_data, bins=SEX_MAX-SEX_MIN+1, range=(SEX_MIN, SEX_MAX+1),
             density=True, edgecolor='white')
    plt.xticks(_SEX_TICK_LOCS, _SEX_TICK_LABELS, fontsize=12)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.ylim(0, 1.025)
    plt.ylabel('Probability', fontsize=12)
    plt.title('SEX Marginal PDF', fontsize=16)
    plt.grid(axis='y')

    plt.subplot(1,4,2)
    plt.hist(ethnicity_data, bins=ETHNICITY_MAX-ETHNICITY_MIN+1,
             range=(ETHNICITY_MIN, ETHNICITY_MAX+1), density=True,
             edgecolor='white')
    plt.xticks(_ETH_TICK_LOCS, _ETH_TICK_LABELS, fontsize=12)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.ylim(0, 1.025)
    plt.ylabel('Probability', fontsize=12)
    plt.title('ETHNICITY Marginal PDF', fontsize=16)
    plt.grid(axis='y')

    plt.subplot(1,4,3)
    plt.hist(case_status_data, bins=CASE_STATUS_MAX-CASE_STATUS_MIN+1,
             range=(CASE_STATUS_MIN, CASE_STATUS_MAX+1), density=True, edgecolor='white')
    plt.ylabel('Probability', fontsize=12)
    plt.xticks(_CS_TICK_LOCS, _CS_TICK_LABELS, fontsize=12)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.ylim(0, 1.025)
    plt.title('CASE_STATUS Marginal PDF', fontsize=16)
    plt.grid(axis='y')

    plt.subplot(1,4,4)
    plt.hist(pregnant_data, bins=PREGNANT_MAX-PREGNANT_MIN+1,
             range=(PREGNANT_MIN, PREGNANT_MAX+1), density=True,
             edgecolor='white')
    plt.xticks(_PREG_TICK_LOCS, _PREG_TICK_LABELS, fontsize=12)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.ylim(0, 1.025)
    plt.ylabel('Probability', fontsize=12)
    plt.title('PREGNANT Marginal PDF', fontsize=16)
    plt.grid(axis='y')

    #plt.subplots_adjust(hspace=0.6)
    plt.show()


###############################################################################
def plot_pdf_ecdf(state):

    race_max = state.race_max
    county_max = state.county_max

    ecdf_age = state.ecdf_age
    ecdf_sex = state.ecdf_sex
    ecdf_race = state.ecdf_race
    ecdf_ethnicity = state.ecdf_ethnicity
    ecdf_county = state.ecdf_county
    ecdf_case_status = state.ecdf_case_status
    ecdf_pregnant = state.ecdf_pregnant

    age_data = state.age_data
    sex_data = state.sex_data
    race_data = state.race_data
    ethnicity_data = state.ethnicity_data
    county_data = state.county_data
    case_status_data = state.case_status_data
    pregnant_data = state.pregnant_data

    AGE_MIN  = HL7.AGE_MIN
    AGE_MAX  = HL7.AGE_MAX
    SEX_MIN  = HL7.SEX_MIN
    SEX_MAX  = HL7.SEX_MAX
    RACE_MIN = HL7.RACE_MIN
    ETHNICITY_MIN = HL7.ETHNICITY_MIN
    ETHNICITY_MAX = HL7.ETHNICITY_MAX
    CASE_STATUS_MIN = HL7.CASE_STATUS_MIN
    CASE_STATUS_MAX = HL7.CASE_STATUS_MAX
    COUNTY_MIN = HL7.COUNTY_MIN
    PREGNANT_MIN = HL7.PREGNANT_MIN
    PREGNANT_MAX = HL7.PREGNANT_MAX

    # race_max is determined from the data, not from a file of constants
    RACE_MAX = race_max

    # county_max is determined from the data, not from a file of constants
    COUNTY_MAX = county_max

    # FIG 1: AGE, RACE, COUNTY
    plt.figure(figsize=(30,16))

    # marginal distributions (histograms of each variable)

    # age
    plt.subplot(3,3,1)
    plt.hist(age_data, bins=AGE_MAX-AGE_MIN+1, range=(AGE_MIN, AGE_MAX+1),
             density=True, align='left', edgecolor='white')
    plt.xlabel('Age in Years')
    plt.ylabel('Probability')
    plt.title('Age Marginal PDF')

    # race
    plt.subplot(3,3,2)
    if RACE_MAX > 80:
        # don't use a white edgecolor, bins can be erased for large bin numbers
        plt.hist(race_data, bins=RACE_MAX-RACE_MIN+1, range=(RACE_MIN, RACE_MAX+1),
                 density=True)
    else:
        plt.hist(race_data, bins=RACE_MAX-RACE_MIN+1, range=(RACE_MIN, RACE_MAX+1),
                 density=True, edgecolor='white')
    if RACE_MAX < 20:
        # one tick for every value
        plt.xticks(range(0, RACE_MAX+1, 1))
    elif RACE_MAX < 40:
        # one tick every two values
        plt.xticks(range(0, RACE_MAX+1, 2))
    else:
        # use the default
        pass
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.ylim(0, 1.025)
    plt.xlabel('Remapped Race Code')
    plt.ylabel('Probability')
    plt.title('Race Marginal PDF')

    # county
    plt.subplot(3,3,3)
    if COUNTY_MAX > 80:
        # don't use a white edgecolor, bins can be erased for large bin numbers
        plt.hist(county_data, bins=COUNTY_MAX-COUNTY_MIN+1,
                 range=(COUNTY_MIN, COUNTY_MAX+1), density=True)
    else:
        plt.hist(county_data, bins=COUNTY_MAX-COUNTY_MIN+1,
                 range=(COUNTY_MIN, COUNTY_MAX+1), density=True,
                 edgecolor='white')
    plt.xlabel('Remapped County Code')
    plt.ylabel('Probability')
    plt.title('County Marginal PDF')

    # ECDFs

    plt.subplot(3,3,4)
    x_vals = np.linspace(AGE_MIN, AGE_MAX, num=AGE_MAX-AGE_MIN+1)
    y_vals = [ecdf_age(x) for x in x_vals]
    plt.step(x_vals, y_vals, where='post')
    plt.xlim(AGE_MIN-5, AGE_MAX+5)
    plt.ylim(-0.1,1.1)
    plt.xlabel('Age in Years')
    plt.ylabel('Probability')
    plt.title('Age ECDF')
    plt.grid()

    plt.subplot(3,3,5)
    x_vals = np.linspace(RACE_MIN, RACE_MAX, num=RACE_MAX-RACE_MIN+1)
    y_vals = [ecdf_race(x) for x in x_vals]
    plt.step(x_vals, y_vals, where='post')
    if RACE_MAX < 20:
        # one tick for every value
        plt.xticks(range(0, RACE_MAX+1, 1))
    elif RACE_MAX < 40:
        # one tick every two values
        plt.xticks(range(0, RACE_MAX+1, 2))
    else:
        # use the default
        pass
    plt.ylim(-0.1, 1.1)
    plt.xlabel('Remapped Race Code')
    plt.ylabel('Probability')
    plt.title('Race ECDF')
    plt.grid()

    plt.subplot(3,3,6)
    x_vals = np.linspace(COUNTY_MIN, COUNTY_MAX, num=COUNTY_MAX-COUNTY_MIN+1)
    y_vals = [ecdf_county(x) for x in x_vals]
    plt.step(x_vals, y_vals, where='post')
    plt.xlim(COUNTY_MIN-1, COUNTY_MAX+1)
    plt.ylim(-0.1,1.1)
    plt.xlabel('Remapped County Code')
    plt.ylabel('Probability')
    plt.title('County ECDF')
    plt.grid()

    # inverse ECDFs

    plt.subplot(3,3,7)
    n_empirical = len(age_data)
    x_vals = [i * (1.0/n_empirical) for i in range(0, n_empirical+1)]
    y_vals = np.array([ecdf_age.inv(u) for u in x_vals])
    plt.step(x_vals, y_vals)
    plt.ylim([AGE_MIN-5, AGE_MAX+5])
    plt.xlabel('Probability')
    plt.ylabel('Age in Years')
    plt.title('Age Inverse ECDF')
    plt.grid()

    plt.subplot(3,3,8)
    n_empirical = len(race_data)
    x_vals = [i * (1.0/n_empirical) for i in range(0, n_empirical+1)]
    y_vals = np.array([ecdf_race.inv(u) for u in x_vals])
    plt.step(x_vals, y_vals)
    if RACE_MAX < 20:
        # one tick for every value
        plt.yticks(range(0, RACE_MAX+1, 1))
    elif RACE_MAX < 40:
        # one tick every two values
        plt.yticks(range(0, RACE_MAX+1, 2))
    else:
        # use the default
        pass
    plt.ylim(RACE_MIN-1, RACE_MAX+1)
    plt.xlabel('Probability')
    plt.ylabel('Remapped Race Code')
    plt.title('Race Inverse ECDF')
    plt.grid()

    plt.subplot(3,3,9)
    n_empirical = len(county_data)
    x_vals = [i * (1.0/n_empirical) for i in range(0, n_empirical+1)]
    y_vals = np.array([ecdf_county.inv(u) for u in x_vals])
    plt.step(x_vals, y_vals)
    plt.ylim(COUNTY_MIN-1, COUNTY_MAX+1)
    plt.xlabel('Probability')
    plt.ylabel('Remapped County Code')
    plt.title('County Inverse ECDF')
    plt.grid()

    plt.subplots_adjust(hspace=0.4, wspace=0.25)
    plt.show()
    
    # FIG 2: SEX, ETHNICITY, CASE_STATUS, PREGNANT

    plt.figure(figsize=(30,16))

    plt.subplot(3,4,1)
    plt.hist(sex_data, bins=SEX_MAX-SEX_MIN+1, range=(SEX_MIN, SEX_MAX+1),
             density=True, edgecolor='white')
    plt.xticks(_SEX_TICK_LOCS, _SEX_TICK_LABELS, fontsize=12)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.ylim(0, 1.025)
    plt.ylabel('Probability')
    plt.title('Sex Marginal PDF')

    plt.subplot(3,4,2)
    plt.hist(ethnicity_data, bins=ETHNICITY_MAX-ETHNICITY_MIN+1,
             range=(ETHNICITY_MIN, ETHNICITY_MAX+1), density=True,
             edgecolor='white')
    plt.xticks(_ETH_TICK_LOCS, _ETH_TICK_LABELS, fontsize=12)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.ylim(0, 1.025)
    plt.ylabel('Probability')
    plt.title('Ethnicity Marginal PDF')


    plt.subplot(3,4,3)
    plt.hist(case_status_data, bins=CASE_STATUS_MAX-CASE_STATUS_MIN+1,
             range=(CASE_STATUS_MIN, CASE_STATUS_MAX+1), density=True,
             edgecolor='white')
    plt.xticks(_CS_TICK_LOCS, _CS_TICK_LABELS, fontsize=12)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.ylim(0, 1.025)
    plt.ylabel('Probability')
    plt.title('Case Status Marginal PDF')

    plt.subplot(3,4,4)
    plt.hist(pregnant_data, bins=PREGNANT_MAX-PREGNANT_MIN+1,
             range=(PREGNANT_MIN, PREGNANT_MAX+1), density=True,
             edgecolor='white')
    plt.xticks(_PREG_TICK_LOCS, _PREG_TICK_LABELS, fontsize=12)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.ylim(0, 1.025)
    plt.ylabel('Probability', fontsize=12)
    plt.title('Pregnant Marginal PDF', fontsize=16)

    # ECDFs

    plt.subplot(3,4,5)
    x_vals = _SEX_TICK_LOCS
    y_vals = [ecdf_sex(x) for x in x_vals]
    plt.step(x_vals, y_vals, where='post')
    plt.xticks(_SEX_TICK_LOCS, _SEX_ECDF_TICK_LABELS, fontsize=12)
    plt.ylim(-0.1, 1.1)
    plt.ylabel('Probability')
    plt.title('Sex ECDF')
    plt.grid()

    plt.subplot(3,4,6)
    x_vals = _ETH_TICK_LOCS
    y_vals = [ecdf_ethnicity(x) for x in x_vals]
    plt.step(x_vals, y_vals, where='post')
    plt.xticks(_ETH_TICK_LOCS, _ETH_ECDF_TICK_LABELS)
    plt.ylim(-0.1, 1.1)
    plt.ylabel('Probability')
    plt.title('Ethnicity ECDF')
    plt.grid()

    plt.subplot(3,4,7)
    x_vals = _CS_TICK_LOCS
    y_vals = [ecdf_case_status(x) for x in x_vals]
    plt.step(x_vals, y_vals, where='post')
    plt.xticks(_CS_TICK_LOCS, _CS_ECDF_TICK_LABELS)
    plt.ylim(-0.1, 1.1)
    plt.ylabel('Probability')
    plt.title('Case_Status ECDF')
    plt.grid()

    plt.subplot(3,4,8)
    x_vals = _PREG_TICK_LOCS
    y_vals = [ecdf_pregnant(x) for x in x_vals]
    plt.step(x_vals, y_vals, where='post')
    plt.xticks(_PREG_TICK_LOCS, _PREG_TICK_LABELS)
    plt.ylim(-0.1, 1.1)
    plt.ylabel('Probability')
    plt.title('Pregnant ECDF')
    plt.grid()

    # inverse ECDFs

    plt.subplot(3,4,9)
    n_empirical = len(sex_data)
    x_vals = [i * (1.0/n_empirical) for i in range(0, n_empirical+1)]
    y_vals = np.array([ecdf_sex.inv(u) for u in x_vals])
    plt.step(x_vals, y_vals)
    plt.yticks(_SEX_TICK_LOCS, _SEX_ECDF_TICK_LABELS)
    plt.xlabel('Probability')
    plt.title('Sex Inverse ECDF')
    plt.grid()

    plt.subplot(3,4,10)
    n_empirical = len(ethnicity_data)
    x_vals = [i * (1.0/n_empirical) for i in range(0, n_empirical+1)]
    y_vals = np.array([ecdf_ethnicity.inv(u) for u in x_vals])
    plt.step(x_vals, y_vals)
    plt.yticks(_ETH_TICK_LOCS, _ETH_ECDF_TICK_LABELS)
    plt.xlabel('Probability')
    plt.title('Ethnicity Inverse ECDF')
    plt.grid()

    plt.subplot(3,4,11)
    n_empirical = len(case_status_data)
    x_vals = [i * (1.0/n_empirical) for i in range(0, n_empirical+1)]
    y_vals = np.array([ecdf_case_status.inv(u) for u in x_vals])
    plt.step(x_vals, y_vals)
    plt.yticks(_CS_TICK_LOCS, _CS_ECDF_TICK_LABELS)
    plt.xlabel('Probability')
    plt.title('Case_Status Inverse ECDF')
    plt.grid()

    plt.subplot(3,4,12)
    n_empirical = len(pregnant_data)
    x_vals = [i * (1.0/n_empirical) for i in range(0, n_empirical+1)]
    y_vals = np.array([ecdf_pregnant.inv(u) for u in x_vals])
    plt.step(x_vals, y_vals)
    plt.yticks(_PREG_TICK_LOCS, _PREG_ECDF_TICK_LABELS)
    plt.xlabel('Probability')
    plt.title('Pregnant Inverse ECDF')
    plt.grid()

    plt.subplots_adjust(hspace=0.4, wspace=0.25)
    plt.show()


###############################################################################
def plot_tuple_pdf_ecdf(tuple_count,
                        tuple_samples,
                        tuple_ecdf,
                        distributions_only = False):
    """
    Plot the date tuple PDF, ECDF, and inverse ECDF.
    """

    if distributions_only:
        plt.figure(figsize=(30,10))
        num_plots = 1
    else:
        plt.figure(figsize=(30,30))
        num_plots = 3

    plt.subplot(num_plots,1,1)
    plt.hist(tuple_samples, bins=500, density=True, align='left', edgecolor='white')
    plt.xlabel('Date Tuple Index', fontsize=16)
    plt.ylabel('Tuple Probability', fontsize=16)
    plt.title('Date Tuple PDF', fontsize=20)
    plt.grid()

    if not distributions_only:
        plt.subplot(3,1,2)
        x_vals = np.linspace(0, tuple_count, num=1000)
        y_vals = [tuple_ecdf(x) for x in x_vals]
        plt.step(x_vals, y_vals, where='post')
        plt.ylim(-0.1, 1.1)
        plt.xlabel('Tuple Index', fontsize=16)
        plt.ylabel('Probability', fontsize=16)
        plt.title('Date Tuple ECDF', fontsize=20)
        plt.grid()

        plt.subplot(3,1,3)
        x_vals = [i*0.001 for i in range(0, 1001)]
        y_vals = np.array([tuple_ecdf.inv(x) for x in x_vals])
        plt.step(x_vals, y_vals)
        plt.xlabel('Probability', fontsize=16)
        plt.ylabel('Date Tuple Index', fontsize=16)
        plt.title('Date Tuple Inverse ECDF', fontsize=20)
        plt.grid()

    plt.show()


###############################################################################
def get_corr_matrix_elements(r, c, dim, result_dict):
    """
    extract the desired matrix element from all matrices;
    return the truth value and a list of values of that element
    """
    
    results = []
    
    if dim not in result_dict:
        return None, []

    # check to see if the matrix element exists in the dim x dim synthetic tau matrix
    if r >= dim or c >= dim:
        return None, []

    # dict of results for this dimension
    dim_dict = result_dict[dim]

    # data-derived Kendall's tau matrix
    tau_original = dim_dict['tau_original']

    # get truth value from largest matrix for this run
    truth_value = tau_original[r,c]

    # list of run dicts for this dimension
    run_dict_list = dim_dict['runs']
    for rd in run_dict_list:
        # number of samples
        # n = rd['n']
        # synthetic Kendall's tau matrix
        tau_synthetic = rd['tau_synthetic']
        # extract matrix element
        matrix_elt = tau_synthetic[r,c]
        results.append(matrix_elt)
         
    return truth_value, results


###############################################################################
def plot_matrix_elements(variable_names,
                         result_dict,
                         min_dim,
                         max_dim):

    # plot nontrivial Kendall's tau matrix elements vs. number of samples

    ax_xmin = 16
    ax_xmax = 524288
    ax_ymin = -1.0
    ax_ymax = 1.0

    # get n values from the largest matrix example (same for all, though)
    n_vals = []
    dim_result = result_dict[max_dim]
    run_dict_list = dim_result['runs']
    for rd in run_dict_list:
        n_vals.append(rd['n'])

    plt.figure(figsize=(36, 18))

    # size of plot matrix for given max_dim
    plot_rows = {2:1, 3:1, 4:2, 5:2, 6:3}
    plot_cols = {2:1, 3:3, 4:3, 5:5, 6:5}

    plot_index = 1
    for r in range(0, max_dim):
        for c in range(r+1, max_dim):
            # subplot declaration
            plt.subplot(plot_rows[max_dim], plot_cols[max_dim], plot_index)

            # plot the curves
            for d in range(min_dim, max_dim+1):
                # get data for each nval for this tau[r,c] and dimension d
                truth_val, elt_vals = get_corr_matrix_elements(r, c, d, result_dict)
                if len(elt_vals) > 0:
                    plt.plot(n_vals, elt_vals, marker='o', label='{0}D'.format(d))
            plt.hlines(truth_val, ax_xmin, ax_xmax, linestyles='dotted', 
                       color='r', linewidth=1, label='Actual = {0:.3f}'.format(truth_val))
            plt.xscale('log')
            plt.xlabel('Samples')
            plt.ylabel('Element ({0}, {1})'.format(r, c))
            plt.title('{0} and {1}'.format(variable_names[r], variable_names[c]))
            plt.legend()
            plt.axis([ax_xmin, ax_xmax, ax_ymin, ax_ymax])
            plt.grid()
            plot_index += 1
    plt.show()


###############################################################################
def plot_timeseries_result(original, synthetic, zoomed=False):

    max_orig = np.max(original)
    mean_orig = np.mean(original)

    diff = synthetic - original
    max_diff = np.max(diff)

    figsize=(20,24)
    plot_count = 6
    if not zoomed:
        figsize=(20,12)
        plot_count = 3

    plt.figure(figsize=figsize)

    # plot the complete view first
    plt.subplot(plot_count, 1, 1)
    plt.plot(range(len(original)), original, color='green')
    plt.ylabel('COUNT', fontsize=16)
    plt.title('Original', fontsize=16)
    plt.ylim(0, max_orig)
    plt.grid()

    plt.subplot(plot_count, 1, 2)
    plt.plot(range(len(synthetic)), synthetic, color='blue')
    plt.ylabel('COUNT', fontsize=16)
    plt.title('Synthetic', fontsize=16)
    plt.ylim(0, max_orig)
    plt.grid()

    plt.subplot(plot_count, 1, 3)
    plt.plot(range(len(diff)), diff, color='gray')
    plt.hlines(0, 0, len(diff)-1, colors=['red'])
    plt.ylabel('COUNT', fontsize=16)
    plt.title('Diff (Synthetic - Original)', fontsize=16)
    plt.grid()

    if zoomed:
        # zoom in if the timeseries has big spikes
        if max_orig > 10*mean_orig:
            max_orig = 10*mean_orig
        else:
            max_orig = max(max_orig, 10)

        plt.subplot(plot_count, 1, 4)
        plt.plot(range(len(original)), original, color='green')
        plt.ylabel('COUNT', fontsize=16)
        plt.title('Original (Zoomed)', fontsize=16)
        plt.ylim(0, max_orig)
        plt.grid()

        plt.subplot(plot_count, 1, 5)
        plt.plot(range(len(synthetic)), synthetic, color='blue')
        plt.ylabel('COUNT', fontsize=16)
        plt.title('Synthetic (Zoomed)', fontsize=16)
        plt.ylim(0, max_orig)
        plt.grid()

        plt.subplot(plot_count, 1, plot_count)
        plt.plot(range(len(diff)), diff, color='gray')
        plt.hlines(0, 0, len(diff)-1, colors=['red'])
        plt.ylabel('COUNT', fontsize=16)
        plt.title('Diff (Synthetic - Original)', fontsize=16)
        max_y = 0.1 * max_diff
        plt.ylim(-max_y, max_y)
        plt.grid()
    
    plt.show()


###############################################################################
def _align_dates(indices, data_map, min_date):
    
    # extract map data as (datetime, count) tuples
    #data_tuples = [(datetime.datetime.strptime(k[:10], '%Y-%m-%d'), v) for k,v in data_map.items()]
    data_tuples = [(dt, count) for dt,count in data_map.items()]
    # convert to an offset in days from the min date
    data_tuples = [((k-min_date).days, v) for k,v in data_tuples]
    # sort by offset
    data_tuples = sorted(data_tuples, key=lambda x: x[0])

    yvals = [0]*len(indices)
    for k,v in data_tuples:
        # add to plot if in range
        if k >= 0 and k < len(yvals):
            yvals[k] = v
        
    return yvals


###############################################################################
def plot_date_signals(dates, signal, maps):
    """
    """

    # all maps contain datetime objects
    first_elec_submit_dt_map = maps[0]
    diag_dt_map              = maps[1]
    died_dt_map              = maps[2]
    report_dt_new_map        = maps[3]
    hosp_admit_dt_map        = maps[4]
    illness_onset_dt_map     = maps[5]
    invest_start_dt_map      = maps[6]
    
    # offset in days from min_date, include final date
    xvals = range(len(dates))
    assert len(xvals) == len(dates)
    min_date = datetime.datetime.strptime(dates[0], '%Y-%m-%d')

    min_x = 0
    max_x = xvals[-1] + 1

    # plot the signal
    plt.figure(figsize=(20,4))
    plt.plot(xvals, signal)
    plt.xlim([min_x, max_x])
    plt.ylabel('Case Count')
    plt.title('Signal Derived from Reference Date', fontsize=16)
    plt.grid()
    plt.show()

    plt.figure(figsize=(20,20))

    plt.subplot(5, 1, 1)
    yvals = _align_dates(xvals, first_elec_submit_dt_map, min_date)
    plt.plot(xvals, yvals)
    plt.xlim([min_x, max_x])
    plt.ylabel('Case Count')
    plt.title('first_elec_submit_dt', fontsize=16)
    plt.grid()
    plt.subplot(5, 1, 2)
    yvals = _align_dates(xvals, report_dt_new_map, min_date)
    plt.plot(xvals, yvals)
    plt.xlim([min_x, max_x])
    plt.ylabel('Case Count')
    plt.title('report_dt_new', fontsize=16)
    plt.grid()
    plt.subplot(5, 1, 3)
    yvals = _align_dates(xvals, illness_onset_dt_map, min_date)
    plt.plot(xvals, yvals)
    plt.xlim([min_x, max_x])
    plt.ylabel('Case Count')
    plt.title('illness_onset_dt', fontsize=16)
    plt.grid()
    plt.subplot(5, 1, 4)
    yvals = _align_dates(xvals, diag_dt_map, min_date)
    plt.plot(xvals, yvals)
    plt.xlim([min_x, max_x])
    plt.ylabel('Case Count')
    plt.title('diag_dt', fontsize=16)
    plt.grid()
    plt.subplot(5, 1, 5)
    yvals = _align_dates(xvals, invest_start_dt_map, min_date)
    plt.plot(xvals, yvals)
    plt.xlim([min_x, max_x])
    plt.ylabel('Case Count')
    plt.title('invest_start_dt', fontsize=16)
    plt.grid()
    plt.show()
