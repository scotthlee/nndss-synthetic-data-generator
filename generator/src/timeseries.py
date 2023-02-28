#!/usr/bin/env python3
"""
"""

import os
import re
import sys
import math
import cmath
import numpy as np
from scipy import signal as scipy_sig

# set to True to enable debug output
_TRACE = False


###############################################################################
#                                                                             #
#                   Fourier method and related functions                      #
#                                                                             #
###############################################################################


###############################################################################
def enable_debug():
    global _TRACE
    _TRACE = True


###############################################################################
def _to_nonneg_integral(data):
    """
    Convert the 'data' variable (which is a list or a numpy array) to
    nonnegative integer values. Any negative values are replaced by zero,
    and all other values are rounded to the nearest integer. Keep the dtype
    as floating point, since numpy requires floating point values for the
    calculations used by the code below.
    """

    for i,x in enumerate(data):
        if x < 0.0:
            x = 0
        else:
            x = int(x + 0.5)
        data[i] = float(x)


###############################################################################
def _std_dev(samples):
    """
    Incremental algorithm to compute the std dev of a list of samples.

    Adapted from
    https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm.
    """
    
    if len(samples) <= 1:
        return 0.0
    
    m2  = 0
    mean  = 0
    count = 0
    for x in samples:
        count += 1
        delta = x - mean
        mean += delta / count
        delta2 = x - mean
        m2 += delta * delta2

    sample_variance = m2 / (count - 1)
    return math.sqrt(sample_variance)


###############################################################################
def _rolling_std(timeseries, window=7):
    """
    Compute the rolling standard deviation of the given timeseries using a
    window of the given length.
    
    results[0] = std(timeseries[0:w])
    results[1] = std(timeseries[0+1:w+1])
    ...
    
    At the end of the results array, only the samples that fit within the
    window are used.
    """
    
    assert window > 0
    last_index = len(timeseries)-1
    
    results = []
    for i in range(len(timeseries)):
        # don't go past the end of the timeseries
        end = min(i+window, last_index)
        samples = timeseries[i:end]
        results.append(_std_dev(samples))
        
    return results


###############################################################################
def modify_sparse_segments(rng,
                           timeseries,
                           threshold_inc = 0,
                           sparseness_min = 0.9):
    """

    Segment the timeseries and check the sparseness of each segment. If any
    segment is more sparse than the 'SPARSENESS_MIN' threshold value, modify
    some of the nonzero values in the segment. Increase or decrease the value
    by the avg of the nonzero elements in the segment with 50% probability.
    """

    # the segmentation length, in days
    SEGMENT_LEN = 60

    # sparseness threshold (i.e., this fraction of the values in the segment
    # are zero)
    SPARSENESS_MIN = sparseness_min

    # modify the values with this probability
    MODIFICATION_THRESHOLD = 1.0 / 3.0

    ts_len = len(timeseries)
    for i in range(0, ts_len, SEGMENT_LEN):

        # compute sparseness of this segment
        zero_count = 0
        nonzero_sum = 0
        end = min(i+SEGMENT_LEN, ts_len-1)
        if end-i <= 0:
            continue

        for x in timeseries[i:end]:
            if 0 == x:
                zero_count += 1
            else:
                nonzero_sum += x
        sparseness = float(zero_count) / float(end - i)

        if sparseness >= SPARSENESS_MIN:
            nonzero_count = end-i - zero_count
            if nonzero_count < 1:
                continue

            # Compute average value of the nonzero elts in the window.
            # Round to nearest and ensure integer value.
            avg_nonzero = nonzero_sum/float(nonzero_count)
            avg_nonzero = float(int(avg_nonzero + 0.5))
            if avg_nonzero < 1.0:
                avg_nonzero = 1.0
            if avg_nonzero > 3.0:
                avg_nonzero = 3.0

            # modify nonzero values in window with probability given
            # by the modification threshold
            modifications = []
            zero_locations = []
            threshold = min(1.0, MODIFICATION_THRESHOLD + threshold_inc)
            for k in range(i, end):
                x = timeseries[k]
                if 0 == x:
                    # found another zero
                    zero_locations.append(k)
                    continue

                test = rng.uniform(0, 1)
                if test < threshold:
                    # modify value with 50% chance to increase or decrease
                    test = rng.uniform(0, 1)
                    if test >= 0.5:
                        x += avg_nonzero
                    else:
                        x -= avg_nonzero
                        if x < 0:
                            x = 0

                    assert x >= 0
                    modifications.append( (k, x) )

            # randomly shuffle the zero locations within this window
            rng.shuffle(zero_locations)

            # now replace some of the zeros with the modifications
            k = 0
            for index_of_mod, value in modifications:
                index_of_zero = zero_locations[k]
                timeseries[index_of_zero] = value
                timeseries[index_of_mod] = 0
                k += 1

    # ensure at least two nonzero values (needed for kernel density estimation
    # bandwidth calculation)

    nonzero_count = 0
    nonzero_values = set()
    nonzero_indices = set()
    for i,x in enumerate(timeseries):
        if 0 != x:
            nonzero_count +=1
            nonzero_values.add(x)
            nonzero_indices.add(i)

    if 0 == nonzero_count:
        # no synthetic data
        if len(timeseries) >= 1:
            timeseries[0] += 1
        if len(timeseries) >= 2:
            timeseries[1] += 1
        rng.shuffle(timeseries)
    elif 1 == nonzero_count:
        for i in range(len(timeseries)):
            if i in nonzero_indices:
                continue
            timeseries[i] += nonzero_values.pop()
            rng.shuffle(timeseries)
            break

    return timeseries


###############################################################################
def _adaptive_noise(rng, std_dev_timeseries, window=7, scale_factor=1.0):
    """
    Returns samples of zero-mean Gaussian noise, where the stddev of the noise
    is determined from a sliding window.
    """
    assert window > 0

    length = len(std_dev_timeseries)
    last_index = length-1

    results = []
    for i in range(0, length, window):
        # don't go past the end of the timeseries
        end = min(i+window, last_index)
        samples = std_dev_timeseries[i:end]
        
        # compute the mean of the set of std_deviation samples
        avg_std_dev = 0.0
        all_equal = True
        for s in samples:
            avg_std_dev += s
            if s != samples[0]:
                all_equal = False
        if len(samples) > 0:
            avg_std_dev /= float(len(samples))

        # Check for the situation of all samples > 0 but all the same value.
        # Make the stddev nonzero in this case.
        if all_equal and len(samples) > 0 and samples[0] > 0:
            avg_std_dev = 2.0
    
        # generate a window-length number of new noise samples
        noise_std_dev = scale_factor * avg_std_dev
        noise = rng.normal(0, noise_std_dev, size=window)

        results.extend(noise)
        
    return results[:length]


###############################################################################
def _add_noise_to_fft(rng,
                      Y,
                      signal,
                      window,
                      pct = 13.0,
                      std_dev_phase_noise = 0.1,
                      to_nonneg_ints=True):
    """
    Add noise to a given percentage of the high-frequency Fourier coefficients.
    The Fourier transform of the 'signal' timeseries is in Y.
    """

    if _TRACE:
        print('Calling _add_noise_to_fft...')
    
    n = len(Y)
    total_counts = np.sum(signal)

    noisy_Y = np.copy(Y)

    # add noise to Fourier components
    if pct > 0.0:

        # modify Fourier components at this index and higher
        index_start = int( (1.0 - 0.01 * pct) * n)

        # phase of each Fourier component
        new_phases = [math.atan2(y.imag, y.real) for y in Y]

        for i in range(index_start, n):
            # randomly change the phase of this component
            phase_delta = rng.normal(0, std_dev_phase_noise)
            new_phases[i] += phase_delta
            # convert to rectangular form and overwrite the Fourier component
            noisy_Y[i] = cmath.rect(abs(Y[i]), new_phases[i])

    # reconstruct using the noise-modified FFT
    reconstructed = np.fft.ifft(noisy_Y).real

    if to_nonneg_ints:
        _to_nonneg_integral(reconstructed)

    #
    # adaptive Gaussian noise generation
    #

    # remove spikes from the timeseries using a median filter
    filtered = scipy_sig.medfilt(reconstructed, kernel_size=3)
    min_filtered = np.min(filtered)
    max_filtered = np.max(filtered)

    if _TRACE:
        print('\t     min_filtered: {0}'.format(min_filtered))
        print('\t     max_filtered: {0}'.format(max_filtered))

    # some timeseries have data once/week or so, and the filtered
    # series can be all zeros
    if 0 == min_filtered and 0 == max_filtered:
        filtered = reconstructed

    # compute rolling std_deviation of the filtered timeseries
    rolling = _rolling_std(filtered, window)

    if _TRACE:
        print('\t      min_rolling: {0}'.format(np.min(rolling)))
        print('\t      max_rolling: {0}'.format(np.max(rolling)))
    
    # Adjust noise stddev scale factor based on total counts. This is purely
    # empirical, determined by doing lots of runs and studying the appearance
    # of the synthetic timeseries. The amplitude noise adjustments need to
    # increase as the total number of samples decreases. This tends to
    # produce synthetic timeseries with MORE counts than the originals, which
    # helps to ensure that SOMETHING changes in very sparse timeseries.
    if total_counts > 100000:
        scale_factor = 0.1
    elif total_counts > 50000:
        scale_factor = 0.2
    elif total_counts > 40000:
        scale_factor = 0.25
    elif total_counts > 30000:
        scale_factor = 0.3
    elif total_counts > 20000:
        scale_factor = 0.35
    elif total_counts > 10000:
        scale_factor = 0.4
    elif total_counts > 5000:
        scale_factor = 0.8
    elif total_counts > 1000:
        scale_factor = 1.0
    elif total_counts > 600:
        scale_factor = 1.3
    else:
        scale_factor = 1.5

    # generate the adaptive noise and add to the signal
    noise = _adaptive_noise(rng, rolling, window, scale_factor)
    reconstructed += noise
    
    return reconstructed


###############################################################################
def gen_synthetic_fourier(rng,
                          timeseries,
                          pct_to_modify = 30, 
                          std_dev_phase_noise = 0.2, 
                          to_nonneg_ints = True,
                          window = 7):
    """
    Generate a synthetic timeseries via imperfect Fourier reconstruction
    with added adaptive noise. The arguments are:
    
    rng: numpy RNG, create externally with rng = np.random.default_rng(seed)
    timeseries: a list of floats
    pct_to_modify: percentage of the high-frequency Fourier components to
                   modify by adding phase noise
    std_dev_phase_noise: self-descriptive, good value seems to be the default
                         of 0.1
    to_nonneg_ints: if True round all values in the synthetic timeseries to
                    the nearest int
    """

    # compute FFT
    Y = np.fft.fft(timeseries)
    
    # add phase noise to 'pct' of the high-frequency components and reconstruct
    # then add adaptive amplitude noise
    synthetic = _add_noise_to_fft(rng,
                                  Y,
                                  timeseries,
                                  window,
                                  pct=pct_to_modify,
                                  std_dev_phase_noise=std_dev_phase_noise,
                                  to_nonneg_ints=to_nonneg_ints)
    
    if to_nonneg_ints:
        # round to nearest nonnegative integer values (but return floats)
        _to_nonneg_integral(synthetic)
    
    return synthetic
