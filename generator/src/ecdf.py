"""
An inplementation of the empirical cumulative distribution function
and its inverse. The implementation in statsmodels is cryptic and does
not provide an inverse.
"""

import numpy as np
from collections import defaultdict

###############################################################################
class EmpiricalCDF():

    def __init__(self, list_or_dict):
        """
        Constructor, argument is either:

            a list (or numpy.ndarray) of samples

            a dict of value:num_samples_with_this_value pairs
        """
        
        if isinstance(list_or_dict, dict):
            self._init_common(list_or_dict)
        else:
            self._init_list(list_or_dict)


    def _init_common(self, counts):

        x = sorted(counts.keys())
        self.x_min = x[0]
        self.x_max = x[-1]

        # compute values of the ECDF at each data point
        #
        #     ECDF(x) == (number of values <= x)/n

        y = []

        sum_of_counts = 0
        for k,v in counts.items():
            sum_of_counts += v

        n = sum_of_counts
        n_inv = 1.0/float(n)

        cum_sum = 0.0
        for x_val in x:
            count_x = counts[x_val]
            cum_sum += float(count_x)
            y_val = cum_sum * n_inv
            y.append(y_val)

        # the final y-value should be very close to 1.0
        assert np.isclose(1.0, y[-1], rtol = 1.0e-5)

        self.x = np.array(x)
        self.y = np.array(y)
        self.n = n
        

    def _init_list(self, samples):
        assert list == type(samples) or np.ndarray == type(samples)

        # get a count of each data value
        count_dict = defaultdict(int)
        for val in samples:
            count_dict[val] += 1

        self._init_common(count_dict)


    def __call__(self, a):
        """
        Return the empirical CDF evaluated at x=a.
        """
        if a < self.x_min:
            return 0.0
        elif a >= self.x_max:
            return 1.0
        else:
            # find the last index i such that x[i-1] <= a < x[i]
            i = np.searchsorted(self.x, a, side='right') - 1
            return self.y[i]
        

    def inv(self, p):
        """
        Return the inverse of the empirical CDF evaluated at p, which is
        the smallest x such that EmpiricalCDF(x) == p.
        """

        # print('ecdf.inv: ')
        # print('\tp: {0}'.format(p))
        # print('\tself.n: {0}'.format(self.n))
        # print('\tlen(self.x): {0}'.format(len(self.x)))
        # print(self.x)

        assert p >= 0 and p <= 1
        
        index = np.searchsorted(self.y, p)
        assert index >= 0 and index <= self.n
        if index < len(self.x): #self.n:
            return self.x[index]
        else:
            return self.x[-1]


###############################################################################
if __name__ == '__main__':

    # tolerance for floating pt comparisons
    EPS = 1.0e-6
    
    # data containing duplicates
    data = [0,5,5,10,10,15,15,20]
    data_dict = {0:1, 5:2, 10:2, 15:2, 20:1}
    
    n = len(data)
    sum_of_counts = 0
    for k,v in data_dict.items():
        sum_of_counts += v
    assert sum_of_counts == n
    
    ecdf = EmpiricalCDF(data)
    ecdf2 = EmpiricalCDF(data_dict)
    
    #
    # test 1: verify the ECDF value at each data point
    #
    
    # data = [0,     5,     5,     10,    10,    15,    15,    20 ]
    y_gold = [1./8., 3./8., 3./8., 5./8., 5./8., 7./8., 7./8., 1.0]
    for i,x in enumerate(data):
        y = ecdf(x)
        assert abs(y - y_gold[i]) < EPS
        y = ecdf2(x)
        assert abs(y - y_gold[i]) < EPS

    # data = [      0           5           10          15          20        ]
    x_vals = [-0.5, 0.0,  2.5,  5.0,  7.5,  10.0, 12.5, 15.0, 17.5, 20.0, 20.5]
    y_gold = [0.,   1./n, 1./n, 3./n, 3./n, 5./n, 5./n, 7./n, 7./n, 1.,   1.  ]

    #
    # test2: verify the ECDF values at each data point, as well as at values
    #        to each side of each data point
    #
    for i,x in enumerate(x_vals):
        y = ecdf(x)
        assert abs(y - y_gold[i]) < EPS
        y = ecdf2(x)
        assert abs(y - y_gold[i]) < EPS        

    #
    # test 3: verify the inverse ECDF values
    #

    # data = [      0           5           10          15          20        ]
    x_vals = [-0.5, 0.0,  2.5,  5.0,  7.5,  10.0, 12.5, 15.0, 17.5, 20.0, 20.5]
    # y    = [0.,   1./n, 1./n, 3./n, 3./n, 5./n, 5./n, 7./n, 7./n, 1.,   1.  ]
    x_gold = [0.0,  0.0,  0.0,  5.0,  5.0,  10.0, 10.0, 15.0, 15.0, 20.0, 20.0]
    for i,x in enumerate(x_vals):
        y = ecdf(x)
        x1 = ecdf.inv(y)
        assert x1 == x_gold[i]
        y = ecdf2(x)
        x1 = ecdf.inv(y)
        assert x1 == x_gold[i]


    
