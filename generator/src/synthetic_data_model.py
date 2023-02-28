import os
import sys
import json
import math
import datetime
import numpy as np

from scipy import stats

from . import timeseries
from . import jurisdictions
from . import correlation_matrix

# Check the python version and ensure that it's at least 3.6. We want numpy
# to use the python 'secrets' module to initialize the random number generator.
# By doing so, the random numbers will be cryptographically secure and will
# help to minimize any chances of synthetic data re-identification.
# The 'secrets' module first appeared in python 3.6. 
_v_major = sys.version_info.major
_v_minor = sys.version_info.minor
if _v_major < 3 or (3 == _v_major and _v_minor < 6):
    raise SystemExit('This software requires python version 3.6 or greater.')


# set to True to enable debug information
_TRACE = False

_HEMORRHAGIC_FEVER_CODES = {
    10660, # yellow fever
    11630, # ebola hemorrhagic fever
    11631, # marburg hemorrhagic fever
    11632, # lassa fever
    11637, # machupo hemorrhagic fever
    11638, # junin hemorrhagic fever
    11639, # sabia-associated hemorrhagic fever
    11640, # crimean-congo hemorrhagic fever
    11644, # lujo virus
    11648, # guanarito hemorrhagic fever
}

_SYPHILIS_CODES = {
    10310, # syphilis, total primary and secondary
    10311, # syphilis, primary
    10312, # syphilis, secondary
    10313, # syphilis, early latent (also early non-primary, non-secondary)
    10314, # syphilis, late latent
    10316, # syphilis, congenital
    10319, # syphilis, late with clinical manifestations including late benign
           #           syphilis and cardiovascular syphilis
    10320, # syphilis, unknown duration or late
}

_SYPHILIS_PRIMARY_AND_SECONDARY = {
    10310, # syphilis, total primary and secondary
    10311, # syphilis, primary
    10312, # syphilis, secondary    
}

# grouped condition codes
_CODE_GROUPS = {
    # syphilis codes handled as a special case below
    10049:[10056],
    10051:[10064],
    10052:[10065],
    10053:[10062],
    10054:[10061],
    10056:[10049],
    10057:[10063],
    10058:[10066],
    10059:[10068],
    10061:[10054],
    10062:[10053],
    10063:[10057],
    10064:[10051],
    10065:[10052],
    10066:[10058],
    10068:[10059],
    10078:[10079],
    10079:[10078],
    10081:[10082],
    10082:[10081],
    10257:[10258],
    10258:[10257],
    10548:[10549,10550],
    10549:[10548,10550],
    10550:[10548,10549],
    10660:[11630,11631,11632,11637,11638,11639,11640,11644,11648],
    10680:[11704,11705],
    11000:[50242,50265],
    11630:[10660,11631,11632,11637,11638,11639,11640,11644,11648],
    11631:[10660,11630,11632,11637,11638,11639,11640,11644,11648],
    11632:[10660,11630,11631,11637,11638,11639,11640,11644,11648],
    11637:[10660,11630,11631,11632,11638,11639,11640,11644,11648],
    11638:[10660,11630,11631,11632,11637,11639,11640,11644,11648],
    11639:[10660,11630,11631,11632,11637,11638,11640,11644,11648],
    11640:[10660,11630,11631,11632,11637,11638,11639,11644,11648],
    11644:[10660,11630,11631,11632,11637,11638,11639,11640,11648],
    11648:[10660,11630,11631,11632,11637,11638,11639,11640,11644],
    11704:[10680,11705],
    11705:[10680,11704],
    11726:[50221,50223],
    11736:[50222,50224],
    50221:[11726,50223],
    50222:[11736,50224],
    50223:[11726,50221],
    50224:[11736,50222],
    50242:[11000,50265],
    50265:[11000,50242],
}

_DEFAULT_OUTPUT_DIR_NETSS = 'synthetic_results_netss'
_DEFAULT_OUTPUT_DIR_HL7   = 'synthetic_results_hl7'


###############################################################################
def enable_debug():
    global _TRACE
    _TRACE = True


###############################################################################
def error_exit(rng, msg):
    """
    A fatal error occurred, so print the random number generator seed to 
    stdout (for debugging reproducibility) and display an error msg.
    """

    print('*** FATAL ERROR ***')
    print('\tRNG Seed: {0}'.format(rng._bit_generator._seed_seq.entropy))
    raise SystemExit(msg)


###############################################################################
def init_rng(rng_seed=None):
    """
    Initialize the random number generator from either a provided integer seed
    or from the system entropy pool (if rng_seed is None).

    If no seed is provided, numpy extracts 128 bits of entropy from the system
    entropy pool and initializes the generator in a cryptographically secure
    manner for python versions >= 3.6.
    """

    rng = np.random.default_rng(rng_seed)

    if _TRACE:
        print('RNG seed: {0}'.format(rng._bit_generator._seed_seq.entropy))

    return rng


###############################################################################
def to_proper_jurisdiction(user_jurisdiction):
    """
    Convert the user-entered jurisdiction into a properly capitalized full
    jurisdiction string to match those in the file tree.
    """

    str_jurisdiction = None

    # check the jurisdiction for an abbreviation and convert to full name
    if user_jurisdiction.upper() in jurisdictions.ABBREV_MAP:
        str_jurisdiction = jurisdictions.ABBREV_MAP[user_jurisdiction.upper()]
    elif user_jurisdiction.lower() in jurisdictions.LKEY_MAP:
        # get the proper capitalization for the directory in the file tree
        str_jurisdiction = jurisdictions.LKEY_MAP[user_jurisdiction.lower()]

    if str_jurisdiction is None:
        print('\n*** ERROR ***\n'
              'The jurisdiction "{0}" is invalid.'.format(user_jurisdiction))

    return str_jurisdiction


###############################################################################
def get_grouped_codes(code, syphilis_total=False):
    """
    Given a condition code, find all associated grouped codes. If no other
    codes should be grouped, return the code as a single-element list.

    If the code is for a syphilis condition, the Boolean 'syphilis_total'
    specifies whether all codes for syphilis should be grouped together.
    """

    if code in _SYPHILIS_CODES:
        # check if primary or secondary code
        if code in _SYPHILIS_PRIMARY_AND_SECONDARY and not syphilis_total:
            # return smaller primary and secondary group
            code_list = list(_SYPHILIS_PRIMARY_AND_SECONDARY)
        else:
            # user wants all syphilis codes to be grouped
            code_list = list(_SYPHILIS_CODES)
    elif code in _CODE_GROUPS:
        code_list = _CODE_GROUPS[code]
        # add the code itself to the groupe
        code_list.append(code)
    else:
        # no groups for this code, so it is a group by itself
        code_list = [code]

    code_set = set(code_list)
    return sorted(list(code_set))

        
###############################################################################
def build_input_filepaths(netss_dir, user_jurisdiction, code_list):
    """
    Construct fully-qualified paths to all data files. Returns a list of
    strings if successful or None in case of error. If data for the desired
    codes does not exist in a jurisdiction, an empty list is returned.
    """

    str_jurisdiction = to_proper_jurisdiction(user_jurisdiction)
    if str_jurisdiction is None:
        return []
    
    jurisdiction_dir = os.path.join(netss_dir, str_jurisdiction)
    if not os.path.isdir(jurisdiction_dir):
        print('\n*** ERROR ***\n'
              'The path to the jurisdiction directory "{0}" is invalid.'.
              format(jurisdiction_dir))
        return []

    file_list = []
    for code in code_list:
        filename = str(code) + '.csv'
        filepath = os.path.join(jurisdiction_dir, filename)
        if not os.path.isfile(filepath):
            # all codes may not exist in a given jurisdiction
            continue
        file_list.append(filepath)

    return file_list


###############################################################################
def default_output_dir(netss_outdir=True):
    """
    Returns the default directory into which synthetic results will be written.
    """

    if netss_outdir:
        return _DEFAULT_OUTPUT_DIR_NETSS
    else:
        return _DEFAULT_OUTPUT_DIR_HL7


###############################################################################
def default_output_file_name(user_jurisdiction,
                             code_list,
                             syphilis_total):
    """
    Return the default output file name for the given jurisdiction and
    condition code. The default file format is CSV.
    """

    str_jurisdiction = to_proper_jurisdiction(user_jurisdiction)
    if str_jurisdiction is None:
        return None

    # this jurisdiction should be known to exist by this point
    assert str_jurisdiction in jurisdictions.INV_ABBREV_MAP

    # convert to lowercase abbreviation
    abbrev = jurisdictions.INV_ABBREV_MAP[str_jurisdiction].lower()

    if 1 == len(code_list):
        name = 'synthetic_{0}_{1}.csv'.format(code_list[0], abbrev)
    else:
        name_list = ['synthetic']
        sorted_codes = [str(code) for code in sorted(code_list)]
        
        # check for all syphilis or all hemorrhagic fever codes
        all_syphilis = True
        all_hemorrhagic = True
        for code in sorted_codes:
            int_code = int(code)
            if int_code not in _SYPHILIS_CODES:
                all_syphilis = False
            if int_code not in _HEMORRHAGIC_FEVER_CODES:
                all_hemorrhagic = False

            if not all_syphilis and not all_hemorrhagic:
                break

        if all_syphilis and syphilis_total:
            # shorten the name instead of concatenating all the syphilis codes
            name_list.append('syphilis_total')
        elif all_hemorrhagic:
            # ditto for the hemorrhagic fever codes
            name_list.append('hemorrhagic_total')
        else:
            # all codes are used in the file name
            name_list.extend(sorted_codes)

        name_list.append(abbrev)
        name = '_'.join(name_list)
        name += '.csv'

    return name


###############################################################################
def build_output_filepath(output_dir, file_name):
    """
    Check to see if the output dir exists. If so, construct a fully-qualified
    path to the output file and return it. Return None on error.
    """

    outdir = os.path.abspath(output_dir)
    if not os.path.isdir(outdir):
        print('\n*** ERROR ***\n'
              'The output directory "{0}" does not exist.'.format(outdir))
        return None

    output_file_path = os.path.join(outdir, file_name)
    return output_file_path


###############################################################################
def _cholesky(A):
    """
    Compute the Cholesky factorization of matrix A and return it. If the 
    factorization does not exist return None.
    """

    try:
        L = np.linalg.cholesky(A)
    except:
        L = None

    return L


###############################################################################
def _fix_matrix(A,dim):
    """
    Function to shrink covariance matrices that have negative eigenvalues.
    """
    
    # original correlation matrix, eigenvalues, and eigenvectors
    tau_original = (2.0/np.pi)*np.arcsin(A)
    w,V = np.linalg.eig(tau_original)
    
    # check if all eigenvalues are positive
    if (w > 1.0e-10).all():
        # find minimum alpha when moving towards correlation matrix
        for alpha in np.linspace(0,1,1001):
            S = A + (alpha * (tau_original - A))
            w,V = np.linalg.eig(S)
            if (w >= 0).all():
                if min(w) > 1e-10:
                    break
        
        # create new matrix
        S = A + (alpha * (tau_original - A))
        
        # check that the cholesky factorization can be performed
        try:
            np.linalg.cholesky(S)
        except:
            S = None
        
        return S
    
    else:
        # set all negative eigenvalues to a small positive value
        for i in range(len(w)):
            if w[i] <= 0:
                w[i] = 1.0e-4
        
        diag_w = np.diag(w)
        
        # make sure the inverse exists
        try:
            v_inverse = np.linalg.inv(V)
        except:
            v_inverse = None
            
        # create the new correlation matrix
        if v_inverse is not None:
            new_tau_original = np.dot(V, np.dot(diag_w, v_inverse))
        else:
            new_tau_original = None
        
        # find minimum alpha when moving towards correlation matrix
        for alpha in np.linspace(0,1,1001):
            S = A + (alpha * (new_tau_original - A))
            w,V = np.linalg.eig(S)
            if (w > 0).all():
                if min(w) > 1e-10:
                    break
        
        # create new matrix
        S = A + (alpha * (new_tau_original - A))
        
        # check that the cholesky factorization can be performed and an alpha was found
        if alpha != 1:
            try:
                C = np.linalg.cholesky(S)
            except:
                S = None
            return S
        
        # no alpha found
        else:
            # find alpha using identiy matrix
            for alpha in np.linspace(0,1,1001):
                S = A + (alpha * (np.identity(dim) - A))
                w,V = np.linalg.eig(S)
                if (w > 0).all():
                    if min(w) > 1e-10:
                        break
            
            # create new matrix
            S = A + (alpha * (np.identity(dim) - A))
            
            # check if alpha found
            if alpha != 1:
                # check cholesky factorization
                try:
                    C = np.linalg.cholesky(S)
                except:
                    S = None
                return S


###############################################################################
def copula_n_variable_model(num_samples,
                            variable_names,
                            tau_original,
                            ecdf_list,
                            rng):
    """
    """

    ERROR_RETURN = (None, None)
    
    # dimension of this model
    dim = len(variable_names)

    assert len(ecdf_list) == dim
    
    # ESTIMATE the required covariance matrix from the tau matrix;
    # this is exact for 2D only
    cov = np.sin( (np.pi/2.0) * tau_original )

    # This covariance matrix MUST be positive-semidifinite to generate
    # multivariate Gaussian random numbers. If the Cholesky decomposition
    # exists, the matrix is positive-semidefinite.

    # Setting 'allow_singular' to False for a 'fixed' matrix prevents
    # the scikit-learn code from zeroing any eigenvalues that are less
    # than condition number * abs(max_eigenvalue). That condition will
    # likely NOT hold for the modified eigenvalues in the 'fixed' matrix.
    
    L = _cholesky(cov)
    allow_singular = False
    if L is None:
        print('*** Covariance matrix is not positive-semidefinite: ***')
        print()
        print(cov)
        print()
        print('Attempting fix...')
        cov2 = _fix_matrix(cov, dim)
        if cov2 is None:
            print('\tFix failed, using identity for Gaussian RNG covariance.')
            cov = np.identity(dim)
        else:
            # fix might have worked; check if Cholesky factorization exists
            L = _cholesky(cov2)
            if L is None:
                print('*** WARNING ***: copula_n_variable_model: ')
                print('\tCovariance matrix is not positive-semidefinite.')
                print()
                print(cov)
                print()
                print('\tFix failed, using identity for Gaussian RNG covariance.')
                cov = np.identity(dim)
            else:
                # fix worked, Cholesky factorization exists
                print('\tfix worked')
                cov = cov2
                allow_singular = True

    # generate the zero-mean Gaussian-distributed samples

    # the mean is the zero vector in 'dim' dimensions
    mean = np.zeros(dim)

    # create 'frozen' distributions with the specified mean, cov, and 
    # allow_singular flag
    try:
        mv_normal = stats.multivariate_normal(mean=mean,
                                              cov=cov,
                                              allow_singular=allow_singular)
    except np.linalg.LinAlgError:
        # Try again but allow singular.
        mv_normal = stats.multivariate_normal(mean=mean,
                                              cov=cov,
                                              allow_singular=True)

    # generate the multivariate normal samples
    Z = mv_normal.rvs(size=num_samples, random_state=rng)
    
    # use the one-dimensional Gaussian CDF to convert each component to
    # a uniform distribution on [0,1]

    normal = stats.distributions.norm()

    U_list = []
    for i in range(dim):
        if Z.ndim > 1:
            U_i = [normal.cdf(z) for z in Z[:,i]]
        else:
            U_i = [normal.cdf(z) for z in Z]
        U_list.append(U_i)
        
    # use the inverse empirical CDF for each variable to project to
    # desired marginal distributions; this is the synthetic data
    X_list = []
    for i in range(dim):
        X_i = [ecdf_list[i].inv(u) for u in U_list[i]]
        X_list.append(X_i)

    # compute Kendall's tau matrix for the synthetic data
    sys.stdout.flush()
    print("\tComputing Kendall's tau correlation matrix...")
    sys.stdout.flush()
    tau_synthetic = correlation_matrix.kendalls_tau_matrix(X_list)
                
    return X_list, tau_synthetic


###############################################################################
def to_timeseries(file_data):
    """
    Generate a timeseries (numpy array of COUNT variable values) and return
    to the caller.
    """

    dates = []
    signal = []
    num_samples = 0

    # the date for the counts currently being accumulated
    cur_count = 0
    cur_date = file_data[0].event_date
    dates.append(cur_date)

    for r in file_data:

        # file loader should ensure this
        assert r.count is not None

        num_samples += int(r.count)

        # another of the same date; keep accumulating
        if r.event_date == cur_date:
            cur_count += r.count
        else:
            # finished accumulation for cur_date
            signal.append(cur_count)

            cur_count = r.count
            cur_date = r.event_date
            dates.append(cur_date)

    # append the final count
    signal.append(cur_count)

    if _TRACE:
        print('to_timeseries: ')
        print('\tlen(file_data): {0}'.format(len(file_data)))          
        print('\t   len(signal): {0}'.format(len(signal)))
        print('\t    len(dates): {0}'.format(len(dates)))

    assert len(signal) == len(dates)

    # sum the counts of the signal, ensure identical to num_samples
    checksum = 0
    for s in signal:
        checksum += s
    assert checksum == num_samples

    return num_samples, np.array(signal), dates

