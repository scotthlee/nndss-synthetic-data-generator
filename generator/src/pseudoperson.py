#!/usr/bin/env python3
"""
The HL7 generator makes use of the pseudoperson concept for the generation of
correlated date tuples. A pseudoperson is one of the unique combinations of
values of the SEX and CASE_STATUS variables. Since SEX={U, M, F} and
CASE_STATUS = {U, S, P, C}, the pseudoperson categories are 
{UU, US, UP, UC, MU, MS, MP, ...}

Sometimes the synthetic values for SEX and CASE_STATUS generate a pseudoperson
category that is not present in the original dataset. This function looks for
such instances and attempts to make corrections.
"""

import os
import sys

from src import hl7 as HL7
from src import model_data_hl7 as data
from src import synthetic_data_model as model

_MAX_ATTEMPTS = 5


###############################################################################
def correct_invalid_pseudopersons(sample_count_synthetic,
                                  variable_names,
                                  synthetic_data,
                                  tau_original,
                                  cdf_list,
                                  rng):
    """
    Find pseudopersons not present in the original data and attempt to correct
    them.

    This function also checks for any pregnant males and corrects them as well.

    variable_names: list of HL7.NAME_* strings
    sample_count_synthetic: number of synthetic data samples
    synthetic_data: list of lists, one list for each variable in variable_names
    """
    
    sex_data = None
    sex_index = None
    cs_data  = None
    cs_index = None
    preg_data = None
    preg_index = None
    
    # Find the indices for the SEX, CASE_STATUS, and PREGNANT variables.
    # The combination of the first two determines a pseudoperson.
    for i,name in enumerate(variable_names):
        if HL7.NAME_SEX == name:
            sex_data = synthetic_data[i]
            sex_index = i
        elif HL7.NAME_CASE_STATUS == name:
            cs_data = synthetic_data[i]
            cs_index = i
        elif HL7.NAME_PREGNANT == name:
            preg_data = synthetic_data[i]
            preg_index = i


    invalid_index_set = set()

    if sex_data is not None and cs_data is not None:
        for i in range(sample_count_synthetic):
            if not data.is_valid_pseudoperson(sex_data[i], cs_data[i]):
                invalid_index_set.add(i)

    if sex_data is not None and preg_data is not None:
        for i in range(sample_count_synthetic):
            if HL7.SEX_VALUE_MALE == sex_data[i] and HL7.PREGNANT_VALUE_YES == preg_data[i]:
                invalid_index_set.add(i)

    invalid_indices = sorted(list(invalid_index_set))
    invalid_count = len(invalid_indices)

    invalid_frac = float(invalid_count) / float(sample_count_synthetic)
    print('\tFound {0} invalid pseudopersons, {1:.3f}% of the total'.
          format(invalid_count, 100.0*invalid_frac))
    
    # early exit if no corrections needed
    if 0 == invalid_count:
        return synthetic_data

    # make _MAX_ATTEMPTS attempts at correcting invalid pseudopersons
    for attempts in range(_MAX_ATTEMPTS):        
        # generate more pseudopersons, at least 10x desired number
        # min of 1000, max of sample_count_synthetic
        extra_count = 10 * invalid_count
        if extra_count < 1000:
            extra_count = 1000
        if sample_count_synthetic > 1000 and extra_count > sample_count_synthetic:
            extra_count = sample_count_synthetic

        # run the copula model again, but do NOT reseed the RNG
        extra_synthetic_data, tau_extra = model.copula_n_variable_model(extra_count,
                                                                        variable_names,
                                                                        tau_original,
                                                                        cdf_list,
                                                                        rng)
        if extra_synthetic_data is None:
            # something went wrong; the copula model code will provide info
            raise SystemExit('FATAL ERROR: copula_n_variable_model (extra)')

        extra_sex_data  = extra_synthetic_data[sex_index]
        extra_cs_data   = extra_synthetic_data[cs_index]
        if preg_data is not None:
            extra_preg_data = extra_synthetic_data[preg_index]

        correction_indices = []
        for i in range(extra_count):
            sex = extra_sex_data[i]
            is_valid_pseudo = data.is_valid_pseudoperson(sex, extra_cs_data[i])
            is_nonpreg_male = HL7.SEX_VALUE_MALE != sex or \
                              (preg_data is not None and HL7.SEX_VALUE_MALE == sex and \
                               HL7.PREGNANT_VALUE_YES != preg_data[i])
            if is_valid_pseudo and is_nonpreg_male:
                correction_indices.append(i)
                if len(correction_indices) == invalid_count:
                    # found all the corrections we need
                    break

        num_corrections = len(correction_indices)
        print('\tFound {0} of {1} needed corrections.'.format(num_corrections, invalid_count))

        # replace original data with corrections
        for k in range(len(synthetic_data)):
            # k indexes the variables AGE, SEX, ...
            for i, invalid_index in enumerate(invalid_indices):
                # apply the ith correction
                synthetic_data[k][invalid_index] = extra_synthetic_data[k][correction_indices[i]]

        # have completed the first 'num_corrections' of the total
        invalid_count -= num_corrections
        invalid_indices = invalid_indices[num_corrections:]
        
        if invalid_count <= 0:
            break
            
    return synthetic_data

