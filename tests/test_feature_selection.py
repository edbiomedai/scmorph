import random

import numpy as np
import pytest

from scmorph.pp.correlation import xim

"""
XIM correlation coefficient tests adopted from https://github.com/czbiohub/xicor/blob/master/tests/xi_test.py
From Wikipedia:
Anscombe's quartet comprises four data sets that have nearly
identical simple descriptive statistics, yet have very different distributions
and appear very different when graphed. Each dataset consists of eleven
(x,y) points. They were constructed in 1973 by the
statistician Francis Anscombe to demonstrate both the importance of graphing
data before analyzing it and the effect of outliers and other influential
observations on statistical properties.

Note that the tests for XIM do not include ties, because for small datasets they make
the resulting statistic unstable.
"""


@pytest.fixture
def anscombes_quartet():
    anscombes_quartet = {
        "x_1": [10, 8, 13, 9, 11, 14, 6, 4, 12, 7, 5],
        "x_2": [10, 8, 13, 9, 11, 14, 6, 4, 12, 7, 5],
        "y_1": [8.04, 6.95, 7.58, 8.81, 8.33, 9.96, 7.24, 4.26, 10.84, 4.82, 5.68],
        "y_2": [9.14, 8.14, 8.74, 8.77, 9.26, 8.1, 6.13, 3.1, 9.13, 7.26, 4.74],
    }
    anscombes_quartet = {k: np.array(v) for k, v in anscombes_quartet.items()}
    return anscombes_quartet


@pytest.fixture
def anscombes_xis(anscombes_quartet):
    random.seed(2020)
    np.random.seed(2020)
    xis = {
        "xi_1": xim(anscombes_quartet["x_1"], anscombes_quartet["y_1"], 1)[0, 1],
        "xi_2": xim(anscombes_quartet["x_2"], anscombes_quartet["y_2"], 1)[0, 1],
    }
    return xis


def test_xim_correlations(anscombes_xis):
    assert anscombes_xis["xi_1"] == 0.43478260869565233
    assert anscombes_xis["xi_2"] == 0.6086956521739131
