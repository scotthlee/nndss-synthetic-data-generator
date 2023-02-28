#!/usr/bin/env python3
"""
Command-line interface for NETSS synthetic data generation.

To show command-line help:

    python3 ./gen_synthetic_data_netss.py --help

Detailed explanations for each argument can be found in the Jupyter notebook.

REQUIRED ARGUMENTS:

    --netss_dir:    path to the root of the NETSS file tree
                    surround in double quotes for bash shell

    --jurisdiction: string, either an abbreviation or full name
                    surround in double quotes for bash shell

    --code:         five-digit integer condition code
  
OPTIONAL ARGUMENTS:

    --outfile:      fully-qualified path to the output file (including extension);
                    supported extensions are .csv and .json;
                    surround in double quotes for bash shell;
                    if omitted, a default name will be generated

    --num_samples:  integral number of synthetic samples to generate
                    if omitted, the number of samples needed to span the
                    date range of the original data will be generated

    --rng_seed:     integer seed value for the random number generator
                    if omitted, the RNG will be seeded from the system time

FLAGS (these take no arguments):

    --debug            whether to enable debug output
    --disable_grouping whether to disable datafile grouping or not
    --syphilis_total   (syphilis codes only) whether to group all Syphilis
                       datasets together and generate a synthetic dataset
                       from the combined data


EXAMPLES (args shown on separate lines for clarity):

   *** NO SPACES ON EITHER SIDE OF THE EQUAL SIGNS ***


1. Lyme disease (code 11080) in Connecticut, with results written to the
   system-generated filename "synthetic_results_netss/synthetic_11080_Connecticut.csv":

   python3 ./gen_synthetic_data_netss.py
       --netss_dir="/data/csels/preprocessed_data/netss"
       --jurisdiction="CT"
       --code=11080


2. Lyme disease in Connecticut with results written to a specified output file:

   python3 ./gen_synthetic_data_netss.py
       --netss_dir="/data/csels/preprocessed_data/netss"
       --jurisdiction="CT"
       --code=11080
       --outfile="lyme_disease_ct.csv"


3. Syphilis in Oregon, group all syphilis files together,
   write synthetic data to "syphilis_grouped_or.csv":

   python3 ./gen_synthetic_data_netss.py
       --netss_dir="/data/csels/preprocessed_data/netss"
       --jurisdiction="Oregon"
       --code=10311
       --outfile="syphilis_grouped_or.csv"
       --syphilis_total


4. Syphilis in Oregon, group all syphilis files together, 65536 samples,
   write synthetic data to "syphilis_grouped_or_65k.json":

   python3 ./gen_synthetic_data_netss.py
       --netss_dir="/data/csels/preprocessed_data/netss"
       --jurisdiction="Oregon"
       --code=10311
       --outfile="syphilis_grouped_or_65k.json"
       --syphilis_total
       --num_samples=65536


"""

import os
import sys
import time
import argparse
import numpy as np

# support modules
from src import plots
from src import timeseries
from src import jurisdictions
from src import correlation_matrix
from src import model_data_netss as data
from src import synthetic_data_model as model

_VERSION_MAJOR = 0
_VERSION_MINOR = 7

_EXT_CSV  = '.csv'
_EXT_JSON = '.json'


###############################################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Command-line interface for synthetic data generation.')
    
    parser.add_argument('--netss_dir',
                        dest='netss_dir',
                        required=True,
                        help='[REQUIRED] path to the root of the NETSS file tree')
    parser.add_argument('--jurisdiction',
                        dest='jurisdiction',
                        required=True,
                        help='[REQUIRED] jurisdiction of the source data')
    parser.add_argument('--code',
                        dest='condition_code',
                        required=True,
                        help='[REQUIRED] five-digit integer condition code')
    parser.add_argument('--outfile',
                        dest='outfile',
                        help='output filepath')
    parser.add_argument('--num_samples',
                        dest='num_samples',
                        type=int,
                        help='number of samples to generate')
    parser.add_argument('--rng_seed',
                        dest='rng_seed',
                        type=int,
                        help='seed for the random number generator')
    parser.add_argument('--debug',
                        action='store_true',
                        help='print debugging information to stdout')
    parser.add_argument('--disable_grouping',
                        dest='disable_grouping',
                        action='store_true',
                        help='disable source data grouping')
    parser.add_argument('--syphilis_total',
                        dest='syphilis_total',
                        action='store_true',
                        help='combine all syphilis data into a single dataset')

    args = parser.parse_args()

    # get the required args
    netss_dir      = args.netss_dir
    jurisdiction   = args.jurisdiction
    condition_code = args.condition_code

    # set defaults for optional args and flags
    outfile          = None
    num_samples      = None
    rng_seed         = None
    debug            = False
    disable_grouping = False
    syphilis_total   = False

    if args.rng_seed:
        rng_seed = args.rng_seed

    if args.disable_grouping:
        disable_grouping = True

    if args.syphilis_total:
        syphilis_total = True

    if args.outfile:
        outfile = args.outfile
        if not outfile.endswith('csv') and not outfile.endswith('json'):
            print('\n*** Supported output file extensions are .csv and .json ***')
            sys.exit(-1)

    print()
    print('Command-line arguments: ')
    print('\t NETSS directory: {0}'.format(netss_dir))
    print('\t    Jurisdiction: {0}'.format(jurisdiction))
    print('\t  Condition code: {0}'.format(condition_code))
    print('\t     Output file: {0}'.format(outfile))
    print('\t     Num samples: {0}'.format(num_samples))
    print('\t        RNG seed: {0}'.format(rng_seed))
    print('\tDisable grouping: {0}'.format(disable_grouping))
    print('\t  Syphilis total: {0}'.format(syphilis_total))
    print()

    #
    # prepare to run
    #

    if args.debug:
        model.enable_debug()
        data.enable_debug()
        timeseries.enable_debug()

    # check that the NETSS directory exists
    if not os.path.isdir(netss_dir):
        print('\n*** The specified NETSS data directory "{0}" does not exist. ***\n'.
              format(netss_dir))
        sys.exit(-1)

    # get the code list
    condition_code = int(condition_code)
    code_list = [condition_code]
    if not disable_grouping:
        code_list = model.get_grouped_codes(condition_code, syphilis_total)
    print('Code group: {0}'.format(code_list))

    # construct the full paths to all input files
    infile_list = model.build_input_filepaths(netss_dir, jurisdiction, code_list)
    if 0 == len(infile_list):
        print('*** No data in jurisdiction "{0}" for codes {1}. ***'.
              format(jurisdiction, code_list))
        sys.exit(0)
    print('Input files:')
    for f in infile_list:
        print('\t{0}'.format(f))

    # build the output file path

    output_dir = None
    output_file_name = None
    if outfile is None:
        # user did not specify anything on the command line, so use the
        # default name and default CSV format
        output_dir = model.default_output_dir()
        output_file_name = model.default_output_file_name(jurisdiction,
                                                          code_list,
                                                          syphilis_total)
    else:
        # check the file extension and set to CSV if none was specified
        fullname, ext = os.path.splitext(outfile)
        if ext is None:
            # no extension specified, so default to CSV format
            outfile += _EXT_CSV
        else:
            ext = ext.lower()
            if _EXT_CSV != ext and _EXT_JSON != ext:
                print('*** Unsupported format for output file: "{0}"'.
                      format(ext))
                print('Supported extension are {0} and {1}'.
                      format(_EXT_CSV, _EXT_JSON))
                sys.exit(-1)

        # get the output dir, if any, from what the user provided
        output_dir, output_file_name = os.path.split(outfile)
        if output_dir is None:
            output_dir = model.default_output_dir()

    assert output_dir is not None
    assert output_file_name is not None

    output_file_path = model.build_output_filepath(output_dir,
                                                   output_file_name)
    if output_file_path is None:
        print('*** FATAL ERROR: buld_output_filepath ***')
        sys.exit(-1)
    
    print('Output file: \n\t{0}\n'.format(output_file_path))

    # initialize the RNG
    rng = model.init_rng(rng_seed)
    
    # check the sample count, if any
    if num_samples is not None:
        if num_samples <= 0:
            print('\n*** The --num_samples argument "{0}" is invalid. ***\n'.
                  format(num_samples))
            print('Please specify a positive integer.\n')
            sys.exit(-1)
    
    #
    # start the run
    #

    start_time = time.time()

    # load input files
    variable_names = ['AGE', 'SEX', 'RACE', 'HISPANIC', 'CASSTAT', 'COUNTY']
    variable_names, tau_original, ecdf_list, file_data = data.init_model_data(infile_list,
                                                                              variable_names,
                                                                              rng)
    if variable_names is None or 6 != len(variable_names):
        model.error_exit(rng, 'init_model_data')

    # 
    # signal processing
    # 

    signal_processing_start_time = time.time()

    sample_count_orig, signal, dates = model.to_timeseries(file_data)
    print('Timeseries information: ')
    print('\tLength of timeseries: {0}'.format(len(signal)))
    
    # save the time when the Fourier code started
    fourier_start_time = time.time()

    # generate the Fourier result using the default params
    synthetic_fourier = timeseries.gen_synthetic_fourier(rng, timeseries = signal)

    # ensure the synthetic result spans as many days as the original
    assert len(synthetic_fourier) == len(signal)

    # save the time when the Fourier code finished
    fourier_end_time = time.time()
    fourier_elapsed = fourier_end_time - fourier_start_time

    # perform modifications to very sparse segments
    attempts = 0
    threshold_inc = 0.1
    while attempts < 7:
        delta = attempts * threshold_inc
        synthetic_timeseries = timeseries.modify_sparse_segments(rng, synthetic_fourier, delta)

        # compute difference signal
        diff_signal = synthetic_timeseries - signal

        # try again if no difference
        if not np.any(diff_signal):
            attempts += 1
            continue
        else:
            break

    # this is the end of signal processing
    signal_processing_end_time = time.time()

    signal_processing_elapsed = signal_processing_end_time - signal_processing_start_time
    print('Signal processing elapsed time: {0:.3f}s'.format(signal_processing_elapsed))
    print('\tFourier synthesis time: {0:.3f}s'.format(fourier_elapsed))


    #
    # copula model
    #


    # sum the counts in the synthetic timeseries
    sample_count_synthetic = int(np.sum(synthetic_timeseries))
    print(' Total counts in the original timeseries: {0}'.format(sample_count_orig))
    print('Total counts in the synthetic timeseries: {0}'.format(sample_count_synthetic))

    if num_samples is not None:
        # override with user-specified value
        sample_count_synthetic = num_samples

    print('Running copula model...')

    # save the time when the copula model started
    copula_start_time = time.time()

    # run the copula model
    synthetic_data, tau_synthetic = model.copula_n_variable_model(sample_count_synthetic,
                                                                  variable_names,
                                                                  tau_original,
                                                                  ecdf_list,
                                                                  rng)
    if synthetic_data is None:
        # something went wrong; the copula model code will provide info
        model.error_exit(rng, 'copula_n_variable_model')

    # remap the synthetic catagorical values from consecutive integers
    # to the NETSS values
    synthetic_results = data.remap_synthetic_results(synthetic_data,
                                                     variable_names)

    # save the end time
    copula_end_time = time.time()
    print('\tElapsed time: {0:.3f}s'.
          format(copula_end_time - copula_start_time))


    #
    # write output file
    #

    print('Writing output file "{0}"'.format(output_file_path))
    data.write_output_file(output_file_path,
                           sample_count_synthetic,
                           synthetic_timeseries,
                           dates,
                           variable_names,
                           synthetic_results)
    
    end_time = time.time()    
    print('Total elapsed time: {0:.3f}s'.format(end_time - start_time))
    print()

