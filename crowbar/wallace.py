import collections
import pandas as pd
import numpy as np
from typing import Union

mismatch_matrix = collections.namedtuple('mismatch_matrix',
                                         ['a', 'b', 'c', 'd', 'n'])


def contingency_table(calls: pd.DataFrame,
                      partition_a: str, partition_b: str) -> pd.DataFrame:

    ct = pd.crosstab(calls[partition_a], calls[partition_b])

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
        return sum(classifications == i)

    n = len(classifications)
    s = set(classifications)

    return 1 - (sum(n_i(i) * (n_i(i) - 1) for i in s) / (n * (n - 1)))


def wallace(mismatches) -> float:

    return mismatches.a / (mismatches.a + mismatches.b)


def adj_wallace(partition_a: str, partition_b: str,
                calls: pd.DataFrame) -> float:

    ct = contingency_table(calls, partition_a, partition_b)

    mismatches = mismatch(ct)

    sid_b = simpsons(calls[partition_b])

    wallace_a_b = wallace(mismatches)

    wallace_i = 1 - sid_b

    awc = (wallace_a_b - wallace_i) / (1 - wallace_i)

    return awc

