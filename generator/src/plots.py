import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from src import netss as NETSS


###############################################################################
def plot_marginals(state):

    county_max    = state.county_max
    age_data      = state.age_data
    sex_data      = state.sex_data
    race_data     = state.race_data
    hispanic_data = state.hispanic_data
    county_data   = state.county_data
    casstat_data  = state.casstat_data

    AGE_MIN  = NETSS.AGE_MIN
    AGE_MAX  = NETSS.AGE_MAX
    SEX_MIN  = NETSS.SEX_MIN
    SEX_MAX  = NETSS.SEX_MAX
    RACE_MIN = NETSS.RACE_MIN
    RACE_MAX = NETSS.RACE_MAX
    HISPANIC_MIN = NETSS.HISPANIC_MIN
    HISPANIC_MAX = NETSS.HISPANIC_MAX
    CASSTAT_MIN = NETSS.CASSTAT_MIN
    CASSTAT_MAX = NETSS.CASSTAT_MAX
    COUNTY_MIN = NETSS.COUNTY_MIN

    # county_max is determined from the data, not from a file of constants
    COUNTY_MAX = county_max

    # width, height    
    plt.figure(figsize=(20, 10))

    #plt.subplot(321)
    plt.hist(age_data, bins=AGE_MAX-AGE_MIN+1, range=(AGE_MIN, AGE_MAX+1),
             density=True, align='left', edgecolor='white')
    plt.xlabel('Age in Years', fontsize=12)
    plt.ylabel('Probability', fontsize=12)
    plt.title('AGE Marginal PDF', fontsize=16)
    plt.grid(zorder=-5)

    plt.figure(figsize=(20, 10))

    #plt.subplot(325)
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

    plt.figure(figsize=(20,6))

    plt.subplot(141)
    plt.hist(sex_data, bins=SEX_MAX-SEX_MIN+1, range=(SEX_MIN, SEX_MAX+1),
             density=True, edgecolor='white')
    plt.xticks([0, 0.5, 1, 1.5, 2, 2.5, 3], ('', 'U', '', 'M', '', 'F', ''),
               fontsize=12)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.ylim(0, 1.025)
    plt.ylabel('Probability', fontsize=12)
    plt.title('SEX Marginal PDF', fontsize=16)
    plt.grid(axis='y')

    plt.subplot(142)
    plt.hist(race_data, bins=RACE_MAX-RACE_MIN+1, range=(RACE_MIN, RACE_MAX+1),
             density=True, edgecolor='white')
    plt.xticks([0,  0.5, 1,  1.5,  2,  2.5, 3,  3.5, 4,  4.5, 5,  5.5, 6],
               ('', 'U', '', 'NA', '', 'A', '', 'B', '', 'W', '', 'O', ''),
               fontsize=12)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.ylim(0, 1.025)
    plt.ylabel('Probability', fontsize=12)
    plt.title('RACE Marginal PDF', fontsize=16)
    plt.grid(axis='y')

    plt.subplot(143)
    plt.hist(hispanic_data, bins=HISPANIC_MAX-HISPANIC_MIN+1,
             range=(HISPANIC_MIN, HISPANIC_MAX+1), density=True,
             edgecolor='white')
    plt.xticks([0, 0.5, 1, 1.5, 2, 2.5, 3], ('', 'U', '', 'Y', '', 'N', ''),
               fontsize=12)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.ylim(0, 1.025)
    plt.ylabel('Probability', fontsize=12)
    plt.title('HISPANIC Marginal PDF', fontsize=16)
    plt.grid(axis='y')

    plt.subplot(144)
    plt.hist(casstat_data, bins=CASSTAT_MAX-CASSTAT_MIN+1,
             range=(CASSTAT_MIN, CASSTAT_MAX+1), density=True, edgecolor='white')
    plt.ylabel('Probability', fontsize=12)
    plt.xticks([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4],
               ('', 'U', '', 'C', '', 'P', '', 'S', ''), fontsize=12)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.ylim(0, 1.025)
    plt.title('CASSTAT Marginal PDF', fontsize=16)
    plt.grid(axis='y')

    #plt.subplots_adjust(hspace=0.6)
    plt.show()


###############################################################################
def plot_pdf_ecdf(state):

    county_max    = state.county_max
    ecdf_age      = state.ecdf_age
    ecdf_sex      = state.ecdf_sex
    ecdf_race     = state.ecdf_race
    ecdf_hispanic = state.ecdf_hispanic
    ecdf_county   = state.ecdf_county
    ecdf_casstat  = state.ecdf_casstat
    age_data      = state.age_data
    sex_data      = state.sex_data
    race_data     = state.race_data
    hispanic_data = state.hispanic_data
    county_data   = state.county_data
    casstat_data  = state.casstat_data

    AGE_MIN  = NETSS.AGE_MIN
    AGE_MAX  = NETSS.AGE_MAX
    SEX_MIN  = NETSS.SEX_MIN
    SEX_MAX  = NETSS.SEX_MAX
    RACE_MIN = NETSS.RACE_MIN
    RACE_MAX = NETSS.RACE_MAX
    HISPANIC_MIN = NETSS.HISPANIC_MIN
    HISPANIC_MAX = NETSS.HISPANIC_MAX
    CASSTAT_MIN = NETSS.CASSTAT_MIN
    CASSTAT_MAX = NETSS.CASSTAT_MAX
    COUNTY_MIN = NETSS.COUNTY_MIN

    # county_max is determined from the data, not from a file of constants
    COUNTY_MAX = county_max

    plt.figure(figsize=(30,16))

    # marginal distributions (histograms of each variable)

    plt.subplot(3,6,1)
    plt.hist(age_data, bins=AGE_MAX-AGE_MIN+1, range=(AGE_MIN, AGE_MAX+1),
             density=True, align='left', edgecolor='white')
    plt.xlabel('Age in Years')
    plt.ylabel('Probability')
    plt.title('Age Marginal PDF')

    plt.subplot(3,6,2)
    plt.hist(sex_data, bins=SEX_MAX-SEX_MIN+1, range=(SEX_MIN, SEX_MAX+1),
             density=True, edgecolor='white')
    plt.xticks([0, 0.5, 1, 1.5, 2, 2.5, 3], ('', 'U', '', 'M', '', 'F', ''),
               fontsize=12)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.ylim(0, 1.025)
    plt.ylabel('Probability')
    plt.title('Sex Marginal PDF')

    plt.subplot(3,6,3)
    plt.hist(race_data, bins=RACE_MAX-RACE_MIN+1, range=(RACE_MIN, RACE_MAX+1),
             density=True, edgecolor='white')
    plt.xticks([0,  0.5, 1,  1.5,  2,  2.5, 3,  3.5, 4,  4.5, 5,  5.5, 6],
               ('', 'U', '', 'NA', '', 'A', '', 'B', '', 'W', '', 'O', ''),
               fontsize=12)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.ylim(0, 1.025)
    plt.ylabel('Probability')
    plt.title('Race Marginal PDF')

    plt.subplot(3,6,4)
    plt.hist(hispanic_data, bins=HISPANIC_MAX-HISPANIC_MIN+1,
             range=(HISPANIC_MIN, HISPANIC_MAX+1), density=True,
             edgecolor='white')
    plt.xticks([0, 0.5, 1, 1.5, 2, 2.5, 3], ('', 'U', '', 'Y', '', 'N', ''),
               fontsize=12)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.ylim(0, 1.025)
    plt.ylabel('Probability')
    plt.title('Hispanic Marginal PDF')

    plt.subplot(3,6,5)
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

    plt.subplot(3,6,6)
    plt.hist(casstat_data, bins=CASSTAT_MAX-CASSTAT_MIN+1,
             range=(CASSTAT_MIN, CASSTAT_MAX+1), density=True,
             edgecolor='white')
    plt.xticks([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4],
               ('', 'U', '', 'C', '', 'P', '', 'S', ''), fontsize=12)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.ylim(0, 1.025)
    plt.ylabel('Probability')
    plt.title('CASSTAT Marginal PDF')


    # ECDFs

    plt.subplot(3,6,7)
    x_vals = np.linspace(AGE_MIN, AGE_MAX, num=AGE_MAX-AGE_MIN)
    y_vals = [ecdf_age(x) for x in x_vals]
    plt.step(x_vals, y_vals, where='post')
    plt.xlim(AGE_MIN-5, AGE_MAX+5)
    plt.ylim(-0.1,1.1)
    plt.xlabel('Age in Years')
    plt.ylabel('Probability')
    plt.title('Age ECDF')
    plt.grid()

    plt.subplot(3,6,8)
    x_vals = [0, 0.5, 1, 1.5, 2, 2.5, 3]
    y_vals = [ecdf_sex(x) for x in x_vals]
    plt.step(x_vals, y_vals, where='post')
    plt.xticks([0, 0.5, 1, 1.5, 2, 2.5, 3],
               ('U (0)', '', 'M (1)', '', 'F (2)', '', ''), fontsize=12)
    plt.ylim(-0.1, 1.1)
    plt.ylabel('Probability')
    plt.title('Sex ECDF')
    plt.grid()

    plt.subplot(3,6,9)
    x_vals = [0,  1,  2,  3,  4,  5,  6]
    y_vals = [ecdf_race(x) for x in x_vals]
    plt.step(x_vals, y_vals, where='post')
    plt.xticks([0,   1,    2,   3,   4,   5,  6],
               ('U (0)', 'NA (1)', 'A (2)', 'B (3)', 'W (4)', 'O (5)', ''))
    plt.ylim(-0.1, 1.1)
    plt.ylabel('Probability')
    plt.title('Race ECDF')
    plt.grid()

    plt.subplot(3,6,10)
    x_vals = [0, 0.5, 1, 1.5, 2, 2.5, 3]
    y_vals = [ecdf_hispanic(x) for x in x_vals]
    plt.step(x_vals, y_vals, where='post')
    plt.xticks([0, 0.5, 1, 1.5, 2, 2.5, 3],
               ('U (0)', '', 'M (1)', '', 'F (2)', '', ''))
    plt.ylim(-0.1, 1.1)
    plt.ylabel('Probability')
    plt.title('Hispanic ECDF')
    plt.grid()

    plt.subplot(3,6,11)
    x_vals = np.linspace(COUNTY_MIN, COUNTY_MAX, num=COUNTY_MAX-COUNTY_MIN)
    y_vals = [ecdf_county(x) for x in x_vals]
    plt.step(x_vals, y_vals, where='post')
    plt.xlim(COUNTY_MIN-1, COUNTY_MAX+1)
    plt.ylim(-0.1,1.1)
    plt.xlabel('Remapped County Code')
    plt.ylabel('Probability')
    plt.title('County ECDF')
    plt.grid()

    plt.subplot(3,6,12)
    x_vals = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]
    y_vals = [ecdf_casstat(x) for x in x_vals]
    plt.step(x_vals, y_vals, where='post')
    plt.xticks([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4],
               ('U (0)', '', 'C (1)', '', 'P (2)', '', 'S (3)', '', ''))
    plt.ylim(-0.1, 1.1)
    plt.ylabel('Probability')
    plt.title('CASSTAT ECDF')
    plt.grid()

    # inverse ECDFs

    plt.subplot(3,6,13)
    n_empirical = len(age_data)
    x_vals = [i * (1.0/n_empirical) for i in range(0, n_empirical+1)]
    y_vals = np.array([ecdf_age.inv(u) for u in x_vals])
    plt.step(x_vals, y_vals)
    plt.ylim([AGE_MIN-5, AGE_MAX+5])
    plt.xlabel('Probability')
    plt.ylabel('Age in Years')
    plt.title('Age Inverse ECDF')
    plt.grid()

    plt.subplot(3,6,14)
    n_empirical = len(sex_data)
    x_vals = [i * (1.0/n_empirical) for i in range(0, n_empirical+1)]
    y_vals = np.array([ecdf_sex.inv(u) for u in x_vals])
    plt.step(x_vals, y_vals)
    plt.yticks([0, 0.5, 1, 1.5, 2, 2.5, 3],
               ('U (0)', '', 'M (1)', '', 'F (2)', '', ''))
    plt.xlabel('Probability')
    plt.title('Sex Inverse ECDF')
    plt.grid()

    plt.subplot(3,6,15)
    n_empirical = len(race_data)
    x_vals = [i * (1.0/n_empirical) for i in range(0, n_empirical+1)]
    y_vals = np.array([ecdf_race.inv(u) for u in x_vals])
    plt.step(x_vals, y_vals)
    plt.yticks([0,1,2,3,4,5,6],
               ('U (0)', 'NA (1)', 'A (2)', 'B (3)', 'W (4)', 'O (5)', ''))
    plt.xlabel('Probability')
    plt.title('Race Inverse ECDF')
    plt.grid()

    plt.subplot(3,6,16)
    n_empirical = len(hispanic_data)
    x_vals = [i * (1.0/n_empirical) for i in range(0, n_empirical+1)]
    y_vals = np.array([ecdf_hispanic.inv(u) for u in x_vals])
    plt.step(x_vals, y_vals)
    plt.yticks([0, 0.5, 1, 1.5, 2, 2.5, 3],
               ('U (0)', '', 'Y (1)', '', 'N (2)', '', ''))
    plt.xlabel('Probability')
    plt.title('Hispanic Inverse ECDF')
    plt.grid()

    plt.subplot(3,6,17)
    n_empirical = len(county_data)
    x_vals = [i * (1.0/n_empirical) for i in range(0, n_empirical+1)]
    y_vals = np.array([ecdf_county.inv(u) for u in x_vals])
    plt.step(x_vals, y_vals)
    plt.ylim(COUNTY_MIN-1, COUNTY_MAX+1)
    plt.xlabel('Probability')
    plt.ylabel('Remapped County Code')
    plt.title('County Inverse ECDF')
    plt.grid()

    plt.subplot(3,6,18)
    n_empirical = len(casstat_data)
    x_vals = [i * (1.0/n_empirical) for i in range(0, n_empirical+1)]
    y_vals = np.array([ecdf_casstat.inv(u) for u in x_vals])
    plt.step(x_vals, y_vals)
    plt.yticks([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4],
               ('U (0)', '', 'C (1)', '', 'P (2)', '', 'S (3)', '', ''))
    plt.xlabel('Probability')
    plt.title('CASSTAT Inverse ECDF')
    plt.grid()

    plt.subplots_adjust(hspace=0.4, wspace=0.25)
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
def plot_timeseries_result(original, synthetic):

    max_orig = np.max(original)
    max_orig = max(max_orig, 10)

    diff = synthetic - original

    plt.figure(figsize=(20,12))
    
    plt.subplot(3, 1, 1)
    plt.plot(range(len(original)), original, color='green')
    plt.ylabel('COUNT', fontsize=16)
    plt.title('Original', fontsize=16)
    plt.ylim(0, max_orig)
    plt.grid()

    plt.subplot(3,1,2)
    plt.plot(range(len(synthetic)), synthetic, color='blue')
    plt.ylabel('COUNT', fontsize=16)
    plt.title('Synthetic', fontsize=16)
    plt.ylim(0, max_orig)
    plt.grid()

    plt.subplot(3,1,3)
    plt.plot(range(len(diff)), diff, color='gray')
    plt.hlines(0, 0, len(diff)-1, colors=['red'])
    plt.ylabel('COUNT', fontsize=16)
    plt.title('Diff (Synthetic - Original)', fontsize=16)
    plt.grid()
    
    plt.show()
