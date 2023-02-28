import numpy as np
from scipy import stats


###############################################################################
def kendalls_tau_matrix(data_lists):
    """
    Compute the Kendall's tau correlation matrix for the data in 'data_lists'.
    The argument is a list of lists, all of identical length.
    The lists contain the values for the random variable that each list
    represents.

    If any matrix elements are NaN, those are set to 0. This situation will
    occur whenever all values for a random variable are identical, which causes
    a division by zero, and hence the NaN. This variable is uncorrelated with
    all other variables, hence a zero value is appropriate.

    The code ASSUMES that 'data_lists' is nonempty and that all components
    have identical lengths.

    The call to stats.kendalltau might possibly overflow the p-value
    calculation, which is ignored in the code below.
    """

    dim = len(data_lists)
    tau = np.ones((dim, dim))

    for r in range(dim):
        for c in range(r+1, dim):
            tau_rc, p_value = stats.kendalltau(data_lists[r], data_lists[c])
            if np.isnan(tau_rc):
                tau_rc = 0.0
            tau[r,c] = tau_rc
            tau[c,r] = tau[r,c]     

    return tau

