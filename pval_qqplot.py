import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata


def expected_pval(pval):
    """ Get expected (uniform distributed) p-value.
    Args:
        pval (array-like): observed p-value
    Returns:
        expected p-value
    """
    return rankdata(pval)/len(pval)


def pval_qqplot(y, ax=None, **kwargs):
    """ Plot expected vs. observed p-value.
    args:
        y (array-like) : observed p-value
        ax (matplotlib.Axes): axis plot
    """
    if ax is None:
        fig, ax = plt.subplots()
    x = expected_pval(y)
    y = pd.Series(y)
    if np.count_nonzero(y) != len(y):
        smallest = y[y != 0].sort_values().iat[0]
        y.loc[y[y==0].index] = smallest
    
    x_log, y_log = -np.log10(x), -np.log10(y) 
    ax.scatter(x_log, y_log, **kwargs)
    # draw y=x line
    x_left, x_right = ax.get_xlim()
    ax.plot([0, x_right], [0, x_right], color='black', linestyle='dashed')
    ax.set_xlabel('Expected -log10 p-value')
    ax.set_ylabel('Observed -1og10 p-value')
    return ax
    