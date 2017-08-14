import argparse
import functools
import json
import re
import operator
import subprocess
from io import StringIO
from pathlib import Path
from typing import Callable, Dict, List, Sequence, Set, Tuple, Union
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

# Third-party imports
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline

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

def closest_relatives(query: str, genomes: Sequence[str],
                      distances: np.matrix) -> Set[int]:

    def which(data: Sequence, operator_: Callable,
              compared: NUMERIC) -> List[int]:

        return [i for i, v in enumerate(data) if operator_(v, compared)]

    query_index = genomes.index(query)
    query_distances = distances[query_index]

    non_self_distances = query_distances[:query_index] + \
                         query_distances[query_index + 1:]

    minimum_distance = min(non_self_distances)

    closest_indices = which(query_distances, operator.eq, minimum_distance)

    # If the minimum dist is 0, it will still match to self here,
    # so ensure the query index is not in the return value
    closest = set(closest_indices) - {query_index}

    return closest

def closest_relative_allele(query: str, missing: str, calls: pd.DataFrame,
                            distances: np.matrix) -> List[int]:

    genomes = calls.index

    closest_relative_indices = closest_relatives(query, genomes, distances)

    closest_relative_alleles = calls[missing].iloc[closest_relative_indices]

    return closest_relative_alleles

def exclude_matches(include: Set[int], matches: Sequence[int]) -> List[int]:

    return [match for match in matches if match in include]

def order_on_reference(reference: Path, genes: Path,
                       calls: pd.DataFrame) -> pd.DataFrame:

    gene_locations = []

    for gene in genes.glob('*.fasta'):

        gene_name = gene.stem

        with gene.open('r') as fasta:
            rec = next(SeqIO.parse(fasta, 'fasta'))

            query = str(rec.seq)

            loc = find_locus_location(query, reference)

            gene_locations.append((gene_name, loc))

    ordered_gene_locations = sorted(gene_locations, key=operator.itemgetter(1))

    ordered_gene_names, ordered_locs = zip(*ordered_gene_locations)

    return calls.reindex_axis(ordered_gene_names, axis=1)


def find_locus_location(query: str, subject: Path) -> int:

    blast = ('blastn', '-task', 'megablast', '-subject', str(subject),
             '-outfmt', '10')

    string_result = subprocess.run(blast, universal_newlines=True, check=True,
                                   input=query,
                                   stdout=subprocess.PIPE)

    table_result = pd.DataFrame(StringIO(string_result.stdout))

    # best row should be first
    # sstart, ssend = 8, 9
    start, stop = table_result.iloc[0, [8,9]]

    return min(start, stop)


def count_triplets(gene: str, calls: pd.DataFrame) -> TRIPLET_COUNTS:

    def tree():
        return defaultdict(tree)

    left_col = calls.columns[calls.columns.index(gene) - 1]
    right_col = calls.columns[calls.columns.index(gene) + 1]

    triplet_calls = calls[[left_col, gene, right_col]]

    triplets = tree()

    for idx, row in triplet_calls.iterrows():

        left, centre, right = row

        try:
            triplets[left][right][centre] += 1

        except TypeError:  # if not yet initialized
            triplets[left][right][centre] = 1

    # Convert back to standard dict to avoid any weirdness later
    return json.loads(json.dumps(triplets))

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
