import sys
import functools
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from datetime import datetime
from typing import Dict, Optional
import pandas as pd
import numpy as np


def user_msg(*messages):
    """Wrapper for print() that prints to stderr"""
    print(*messages, file=sys.stderr)


def logtime(name):
    """Function decorator that print to stderr the runtime of the
    decorated function.
    """

    def decorator(func):
        """Interface between wrapper and the outer logtime() function"""

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            """Wraps func and prints to STDERR the runtime of func"""

            msg = 'Elapsed time for {}: {}'

            before = datetime.now()

            result = func(*args, *kwargs)

            after = datetime.now()

            user_msg(msg.format(name, after - before))

            return result
        return wrapper
    return decorator


def row_distance(idx: int, row, calls: pd.DataFrame) -> Dict[int, int]:
    """Returns the Hamming distance of non-missing alleles between two strains.

    Results are returned as a dictionary for distances between the query strain
    (strain1) and all the subject strain (each strain2).

    Called by dist_gene()
    """

    strain1 = row

    def non_missing_hamming(j):
        """Returns the distance between two strains, considering only loci
        which are not missing in either individual.
        """

        strain2 = calls[j]

        return sum([a > 0 and b > 0 and a != b
                    for a, b in zip(strain1, strain2)])

    return {j: non_missing_hamming(j) for j in range(idx + 1, len(calls))}


def hamming_distance_matrix(distance_path: Optional[Path], calls: pd.DataFrame,
                            cores: int) -> np.matrix:

    """Returns a Hamming distance matrix of pairwise strain distances.

    First attempts to load a pre-calculated matrix from `distance_path and
    returns that if possible. If there is no distance matrix located at that
    path, it calculates one and saves it there. If no path is provided,
    calculate the distance matrix and return it without saving to disk.
    """

    def dist_gene() -> np.matrix:
        """Returns a Hamming distance matrix of pairwise strain distances."""

        n_row = len(calls)

        dist_mat = np.matrix([np.zeros(n_row) for _ in range(n_row)], dtype=int)

        calls_mat = calls.as_matrix()

        with ProcessPoolExecutor(max_workers=cores) as ppe:
            futures = {i: ppe.submit(row_distance, i, row, calls_mat)
                       for i, row in enumerate(calls_mat)}

        results = {i: j.result() for i, j in futures.items()}

        for i, js in results.items():

            for j in js:
                dist_mat[i, j] = dist_mat[j, i] = results[i][j]

        return dist_mat

    try:

        distances = pd.read_csv(distance_path,
                                header=0, index_col=0).as_matrix()

    except FileNotFoundError:

        distances = dist_gene()

        pd.DataFrame(distances,
                     index=calls.index,
                     columns=calls.index).to_csv(distance_path)

    except ValueError:

        distances = dist_gene()

    return distances