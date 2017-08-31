# import from one level above
import sys
import os.path

up = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
sys.path.append(up)

import crowbar  # main script

# regular imports
import argparse
import random
from functools import partial
from itertools import chain
from pathlib import Path
from typing import Dict, Tuple
# 3rd Party imports
import pandas as pd
import numpy as np

# Complex types
#RESULTS = Dict[Dict[int, int], List[int], Set[Optional[int]]]
#RECOVERED = Dict[str, Dict[str, RESULTS]]
def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('--cores',
                        type=int,
                        default=1,
                        help='Number of CPU cores to use [1]')

    parser.add_argument('--seed',
                        type=int,
                        default=1,
                        help='Seed to initialize the PRNG [1]')

    parser.add_argument('--truncation-probability',
                        type=float,
                        default=0.0,
                        store='trunc_prob',
                        help='Uniform probability that any given \
                              locus will be truncated [0.0]')

    parser.add_argument('--missing-probability',
                        type=float,
                        default=0.0,
                        store='miss_prob',
                        help='Uniform probability that any given \
                              locus will be rendered missing [0.0]')

    parser.add_argument('--distances',
                        type=Path,
                        required=False,
                        help='Path to pre-calculated distance matrix')

    parser.add_argument('--tempdir',
                        type=Path,
                        default=Path('/tmp'),
                        help='Directory for emphermal working files')

    parser.add_argument('calls',
                        type=Path,
                        help='Table of allele calls')

    parser.add_argument('genes',
                        type=Path,
                        help='Directory of gene multifastas')

    parser.add_argument('jsons',
                        type=Path,
                        help='Directory containing MIST results')

    return parser.parse_args()

def sequential_test():

    # for each calls[strain, gene], truncate or vanish
        # get results
        # return most probable match
    pass

def introduce_random_errors(trunc_prob: float, miss_prob: float,
                            calls: pd.DataFrame, jsondir: Path,
                            seed: int) -> Tuple[pd.DataFrame, Dict[str, str]]:

    random.seed(seed)

    error_calls = calls.copy()

    for strain, row in error_calls.iterrows():

        to_truncate = [random.random() < trunc_prob for _ in row]
        to_vanish = [random.random() < miss_prob for _ in row]

        row[to_truncate] = [-1 for _ in row]
        row[to_vanish] = [0 for _ in row]

        truncations = {gene: truncate(strain, gene, jsondir)
                       for gene in row[to_truncate].index}

    return error_calls, truncations


def truncate(strain: str, gene: str, jsondir: Path, tempdir: Path):

    sequence = crowbar.load_fragment(strain, gene, jsondir)

    # Leave at least a 50-mer on each end
    pivot = random.randint(50, len(sequence) - 50)
    side = random.choice((0, 1))

    halves = sequence[:pivot], sequence[:pivot]

    return halves[side]


def recover_alleles(truncations: Dict[str, str], distances: np.matrix,
                    calls: pd.DataFrame, genes: Path):

    gene_names = calls.columns

    results = {}

    for strain, row in calls.iterrows():
        results[strain] = {}

        # treat truncated and missing differently,
        # as missing loci have no fragment to use
        truncated_genes = row[row == -1].index
        missing_genes = row[row == 0].index

        for gene in chain(truncated_genes, missing_genes):

            ### Nearest Neighbour
            closest_alleles = crowbar.closest_relative_allele(strain, gene,
                                                              calls, distances)


            ### Linkage disequilibrium
            triplets = count_triplets(gene, calls)
            left_col = gene_names[gene_names.index(gene) - 1]
            right_col = gene_names[gene_names.index(gene) + 1]

            left, right = calls[[left_col, right_col]].loc[strain]
            linked = triplets[left][right]



            ### Fragment matching
            try:
                matching_fragments = crowbar.fragment_match(gene,
                                                            truncations[gene],
                                                            genes)

                partial_matches_only = partial(crowbar.exclude_matches,
                                               include=matching_fragments)

                closest_alleles = partial_matches_only(closest_alleles)

                linked = {centre: linked[centre]
                          for centre in partial_matches_only(linked.keys())}
            except KeyError:

                matching_fragments = set()

            results[strain][gene] = {'linkage': linked,
                                     'closest': closest_alleles,
                                     'fragments': matching_fragments}

    return results


def recover(truncation_probability: float, missing_probability: float,
            seed: int, calls: pd.DataFrame, jsondir: Path, genes: Path,
            cores: int):

    distances = crowbar.dist_gene(calls, cores)

    error_calls, truncations = introduce_random_errors(truncation_probability,
                                                       missing_probability,
                                                       calls, jsondir, seed)

    recovered_alleles = recover_alleles(truncations, distances,
                                        error_calls, genes)


def main():

    args = arguments()


if __name__ == '__main__':
    main()
