"""
Kernel density estimation module for the synthetic data generator.
"""

import math
import numpy as np
import scipy.integrate as integrate

from scipy import interpolate
from collections import defaultdict
from statsmodels.nonparametric import bandwidths

from .ecdf import EmpiricalCDF
from .netss import AGE_UNKNOWN_REMAP

# This value is used to transform a random variable X to a new random variable
# Y = log(X + KDE_DELTA). The CSELS unknown age value is -1, so this value must
# be greater than 1 for the logarithm to be a real number.
# The resulting density estimate is very sensitive to delta values near 1, so 
# a better value is further out on the log curve near 7 or so.
KDE_DELTA = 7

_SQRT_TWO_PI = math.sqrt(2.0 * math.pi)


###############################################################################
def _kernel_bandwidth(samples, bw_method='silverman'):
    """
    Compute a kernel bandwidth by one of these methods (provided by the
    statsmodels library): 'silverman', 'scott', 'normal'.
    """
    
    if 'silverman' == bw_method.lower():
        bw = bandwidths.bw_silverman(samples)
    elif 'scott' == bw_method.lower():
        bw = bandwidths.bw_scott(samples)
    elif 'normal' == bw_method.lower():
        bw = bandwidths.bw_normal_reference(samples)
    else:
        raise SystemError('_kernel_bandwidth: unknown bandwidth method "{0}"'.
                          format(bw_method))

    # For small sample counts the bandwidth could be negative. Take corrective
    # action in that case.
    if bw <= 0:
        if len(samples) >= 2:
            # compute the sample std deviation (which could be zero)
            bw = np.std(samples)

    # catchall
    if math.isnan(bw) or bw <= 0:
        bw = 5.0 # five years

    return bw


###############################################################################
def _kde_gaussian(x, n, bins, bw):
    """
    Kernel density estimation using a Gaussian kernel; evaluate at
    a single point x.
    """ 
    
    # precompute a few values, to avoid repeated divisions
    inv_n     = 1.0 / float(n)
    inv_bw    = 1.0 / bw
    inv_denom = 1.0 / (bw * _SQRT_TWO_PI)
        
    # evaluate for each bin and multiply by the weight of that bin,
    # instead of for all samples
    val = 0
    for binval, count in bins.items():
        t = (x - binval) * inv_bw
        num = math.exp(-0.5 * t * t)
        # there are 'count' identical copies of this sample value,
        # which is the bin value
        val += count * num * inv_denom
        
    return val * inv_n


###############################################################################
def bounded_gaussian(original_samples, start, stop, bw_method='silverman'):
    """
    Given samples of a random variable X, perform kernel density estimation
    for X using Gaussian kernels.
    """

    num_orig_samples = len(original_samples)

    # remove unknowns from samples
    samples = np.array([s for s in original_samples if AGE_UNKNOWN_REMAP != s])
    num_samples = len(samples)

    # fraction of unknowns in original sample
    unk_frac = (num_orig_samples - num_samples) / num_orig_samples

    # kde to be evaluated on a linear, evenly-spaced grid between [0, stop]
    num_pts = stop - 0 + 1
    xvals = np.linspace(0, stop, num_pts)

    if unk_frac >= 1.0:
        # all age values are unknown, nothing to compute
        kde = np.zeros(len(xvals))
        return xvals, kde, 1.0
    
    # compute the bandwidth by the specified method
    bw = _kernel_bandwidth(samples, bw_method)
    assert bw > 0
    
    # bin the values and count them (speeds up the evaluation over all x)
    bins = defaultdict(float)
    for s in samples:
        bins[s] += 1.0
    
    # evaluate the kde
    kde = np.array([_kde_gaussian(x, num_samples, bins, bw) for x in xvals])
    
    # compute the area under the curve from [0, stop] and rescale so that
    # the integral value is 1.0-unk_frac
    integral = integrate.simps(kde, xvals)
    if integral is not None and not math.isnan(integral) and integral > 0:
        kde = kde * (1.0 - unk_frac) / integral

    # this constraint now holds: unk_frac + area_under_kde_curve == 1

    return xvals, kde, unk_frac


###############################################################################
def _kde_log(x, num_samples, bins, bw, delta):
    """
    Kernel density estimation using a lognormal kernel;
    evaluate at a single point x.
    The 'samples' array is the raw data.
    """
    
    # precompute to avoid repeated divisions (slow) and multiplications
    inv_n       = 1.0 / float(num_samples)
    inv_bw      = 1.0 / bw
    denom_const = bw * _SQRT_TWO_PI
        
    # evaluate for each bin and multiply by the weight of that bin,
    # instead of for all samples
    val = 0
    for binval, count in bins.items():
        # arguments for the log function must be positive
        assert x + delta > 0.0
        assert binval + delta > 0.0
        t = (math.log(x + delta) - math.log(binval + delta)) * inv_bw
        num = math.exp(-0.5 * t * t)
        denom = ( (x + delta) * denom_const)
        # there are 'count' identical copies of this sample value,
        # which is the bin value
        val += count * num / denom
        
    return val * inv_n


###############################################################################
def bounded_lognormal(original_samples, start, stop, bw_method='silverman',
                      delta=KDE_DELTA):
    """
    Given samples of a random variable X, perform kernel density estimation
    for X via the transformation Y = log(X + delta).
    """

    num_orig_samples = len(original_samples)

    # remove unknowns from samples
    samples = np.array([s for s in original_samples if AGE_UNKNOWN_REMAP != s])
    num_samples = len(samples)
    
    # fraction of unknowns in original sample
    unk_frac = (num_orig_samples - num_samples) / num_orig_samples

    # kde to be evaluated on a linear, evenly-spaced grid between [0, stop]
    num_pts = stop - 0 + 1
    xvals = np.linspace(0, stop, num_pts)

    if unk_frac >= 1.0:
        # all age values are unknown, nothing to compute
        kde = np.zeros(len(xvals))
        return xvals, kde, 1.0

    # estimate the bandwidth for the TRANSFORMED samples:
    #     log_samples = log(samples + delta)
    log_samples = np.log(samples + delta)
    bw = _kernel_bandwidth(log_samples, bw_method)
    assert bw > 0
    
    # bin the samples and count them
    bins = defaultdict(float)
    for s in samples:
        bins[s] += 1.0     
    
    # evaluate the KDE
    kde = np.array([_kde_log(x, num_samples, bins, bw=bw, delta=delta) 
                    for x in xvals])

    # compute the area under the curve from [0, stop] and rescale so that
    # the integral value is 1.0-unk_frac
    integral = integrate.simps(kde, xvals)
    if integral is not None and not math.isnan(integral) and integral > 0:
        kde = kde * (1.0 - unk_frac) / integral

    # this constraint now applies: unk_frac + area_under_kde_curve == 1

    return xvals, kde, unk_frac


###############################################################################
def bounded_hybrid(samples, start, stop, bw_method='silverman',
                   delta=KDE_DELTA):
    """
    """

    xvals1, kde1, unk_frac1 = bounded_lognormal(samples, start, stop, bw_method, delta)
    xvals2, kde2, unk_frac2 = bounded_gaussian (samples, start, stop, bw_method)
    
    # must have identical unknowns in both (default tolerance is 1.0e-5)
    assert np.isclose(unk_frac1, unk_frac2)

    # must also haveevaluated KDE on same grid
    assert np.array_equal(xvals1, xvals2)

    unk_frac = unk_frac1
    xvals = xvals1
    
    kde = []
    for i,x in enumerate(xvals):
        w1 = 100.0-4*x
        if w1 < 0:
            w1 = 0
        w2 = 100.0 - w1
        kde_val = (kde1[i] * w1 + kde2[i] * w2) / (w1 + w2)
        if kde_val < 0:
            kde_val = 0
        kde.append(kde_val)

    kde = np.array(kde)

    # compute the area under the curve from [start, stop] and rescale so
    # that it becomes equal to 1 - unk_frac
    integral = integrate.simps(kde, xvals)
    if integral is not None and not math.isnan(integral) and integral > 0:
        kde = kde * (1.0 - unk_frac) / integral
    
    return xvals, kde, unk_frac


###############################################################################
class KdeECDF():
    """
    """

    def __init__(self, samples, l, r, kernel='lognormal'):
        """
        """

        self.ecdf = self.__ecdf_from_kde(samples, l, r, kernel)


    def __ecdf_from_kde(self, samples, l, r, kernel):
        """
        Kernel options are 'lognormal', 'hybrid', and 'gaussian' if unspecified.
        """

        if 'lognormal' == kernel:
            xvals_kde, yvals_kde, unk_frac = bounded_lognormal(samples, l, r)
        elif 'hybrid' == kernel:
            xvals_kde, yvals_kde, unk_frac = bounded_hybrid(samples, l, r)
        else:
            xvals_kde, yvals_kde, unk_frac = bounded_gaussian(samples, l, r)


        # Use the KDE value to compute simulated sample counts at each xvalue.
        # Use a population of 10000 samples; the actual number does not matter
        # for computing an ECDF.
        SIM_SAMPLES = 10000
        unknown_count = int(unk_frac * SIM_SAMPLES)
        sample_dict = {int(AGE_UNKNOWN_REMAP):unknown_count}
        for i, prob in enumerate(yvals_kde):
            count = int(SIM_SAMPLES * prob)
            x = int(xvals_kde[i])
            sample_dict[x] = count

        return EmpiricalCDF(sample_dict)


    def __call__(self, a):
        return self.ecdf(a)


    def inv(self, p):
        return self.ecdf.inv(p)
