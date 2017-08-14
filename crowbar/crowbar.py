import argparse
import functools
import re
import operator
from pathlib import Path
from typing import Dict, Set, Tuple, Sequence, Callable
from concurrent.futures import ProcessPoolExecutor
# Third-party imports
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
import pandas as pd
import numpy as np

# Complex type constants
TRIPLET_COUNTS = Dict[int, Dict[int, Dict[int, int]]]
NUMERIC = Union[int, float]

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('--cores',
                        type=int,
                        default=1,
                        help='Number of CPU cores to use [1]')

    parser.add_argument('calls',
                        type=Path,
                        help='Table of allele calls')

    parser.add_argument('genes',
                        type=Path,
                        help='Directory of gene multifastas')

    return parser.parse_args()


def row_distance(idx: int, row: Tuple[str, pd.Series],
                 calls: pd.DataFrame) -> Dict[int, int]:

    name, strain1 = row

    return {j: sum(strain1 != calls.iloc[j])
            for j in range(idx + 1, len(calls))}


def dist_gene(calls: pd.DataFrame, cores: int) -> np.matrix:


    n_row = len(calls)

    dist_mat = np.matrix([[0 for _ in range(n_row)] for _ in range(n_row)])

    with ProcessPoolExecutor(max_workers=cores) as ppe:
        futures = {i: ppe.submit(row_distance, i, row, calls)
                   for i, row in enumerate(calls.iterrows())}

    results = {i: j.result() for i, j in futures.items()}

    for i, js in results.items():

        for j in js:
            dist_mat[i, j] = dist_mat[j, i] = results[i][j]

    return dist_mat

def closest_relative(query: str, genes: Sequence[str],
                     distances: np.matrix) -> Set[int]:

    def which(data: Sequence, operator: Callable,
              compared: NUMERIC) -> List[int]:

        return [i for i, v in enumerate(data) if operator(v, compared)]

    query_index = genes.index(query)
    query_distances = distances[query_index]

    non_self_distances = query_distances[:query_index] + \
                         query_distances[query_index + 1:]

    minimum_distance = min(non_self_distances)

    closest_indices = which(query_distances, operator.eq, minimum_distance)

    # If the minimum dist is 0, it will still match to self here,
    # so ensure the query index is not in the return value
    closest = set(closest_indices) - {query_index}

    return closest

def closest_relative_allele():
    pass

def order_on_reference(reference: Path, genes: Path,
                       calls: pd.DataFrame) -> pd.DataFrame:
    pass


def linkage_disequilibrium(gene: str, calls: pd.DataFrame) -> TRIPLET_COUNTS:
    pass


def fragment_match(gene: str, fragment: str, genes: Path) -> Set[int]:

    fragment_pattern = re.compile('^{seq}|{seq}$'.format(seq=fragment))

    glob_pattern = '*{}*.f'.format(gene)
    gene_file, *_ = genes.glob(glob_pattern)

    with gene_file.open('r') as fasta:

        result = set(rec.id for rec in SeqIO.parse(fasta, 'fasta')
                     if re.search(fragment_pattern, rec.seq))

    return result

def main():

    args = arguments()

if __name__ == '__main__':
    main()
