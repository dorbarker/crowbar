import collections
import pandas as pd
import numpy as np
from typing import Union

mismatch_matrix = collections.namedtuple('mismatch_matrix',
                                         ['a', 'b', 'c', 'd', 'n'])

wallace_values = collections.namedtuple('wallace_values', ['a_b', 'b_a'])

import warnings
warnings.simplefilter('error')
np.seterr(all='raise')

def contingency_table( partition_a, partition_b) -> pd.DataFrame:

    ct = pd.crosstab(partition_a, partition_b)

    return ct


def mismatch(ct: pd.DataFrame):

    def f(x: Union[pd.DataFrame, pd.Series]) -> np.float64:
        return (x * (x - 1)).values.sum() / 2

    col_sums = ct.sum(0)
    row_sums = ct.sum(1)

    n = sum(col_sums)

    a = f(ct)

    b = f(row_sums) - a

    c = f(col_sums) - a

    d = ((n * (n - 1)) / 2) - (b + a) - c

    return mismatch_matrix(a, b, c, d, n)


def simpsons(classifications: pd.Series) -> float:

    def n_i(i) -> int:
        return np.equal(classifications, i).sum()

    n = len(classifications)
    s = set(classifications)

    return 1 - (sum(n_i(i) * (n_i(i) - 1) for i in s) / (n * (n - 1)))


def wallace(mismatches) -> float:

    return mismatches.a / (mismatches.a + mismatches.b)


def adj_wallace(partition_a, partition_b) -> float:

    ct = contingency_table(partition_a, partition_b)

    mismatches = mismatch(ct)

    sid_b = simpsons(partition_b)

    wallace_a_b = wallace(mismatches)

    wallace_i = 1 - sid_b

    try:
        awc = (wallace_a_b - wallace_i) / (1 - wallace_i)
    except FloatingPointError:
        awc = 1.0


    return awc

def bidirectional_adj_wallace(partition_a, partition_b):

    awc_a_b = adj_wallace(partition_a, partition_b)

    awc_b_a = adj_wallace(partition_b, partition_a)

    return wallace_values(awc_a_b, awc_b_a)

