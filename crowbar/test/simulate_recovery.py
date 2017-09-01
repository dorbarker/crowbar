# regular imports
import sys
import os.path
import argparse
import json
import random
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path
from typing import Dict, Tuple, Union
# 3rd Party imports
import pandas as pd
import numpy as np

up = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
sys.path.append(up)

import crowbar  # main script

# Complex types
TRUNCATIONS = Dict[str, Dict[str, str]]

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

ERROR_CALLS_TRUNCS = Tuple[pd.DataFrame, TRUNCATIONS]
def random_errors(trunc_prob: float, miss_prob: float, calls: pd.DataFrame,
                  jsondir: Path, seed: int) -> ERROR_CALLS_TRUNCS:

    random.seed(seed)

    error_calls = calls.copy()

    truncations = {}

    for strain, row in error_calls.iterrows():

        to_truncate = [random.random() < trunc_prob for _ in row]
        to_vanish = [random.random() < miss_prob for _ in row]

        row[to_truncate] = [-1 for _ in row]
        row[to_vanish] = [0 for _ in row]

        truncations[strain] = {gene: truncate(strain, gene, jsondir)
                               for gene in row[to_truncate].index}


    return error_calls, truncations


def truncate(strain: str, gene: str, jsondir: Path) -> str:

    def load_fragment(strain: str, gene: str, jsondir: Path) -> str:

        """Loads a partial nucleotide alignment from MIST-generated JSON
        output.

        Currently, this function assumes that there is only one cgMLST test per
        JSON file.
        """

        jsonpath = (jsondir / strain).with_suffix('.json')

        with jsonpath.open('r') as json_obj:
            data = json.load(json_obj)

        test_results = data['Results'][0]['TestResults']

        # presumed to be a single test name per file
        test_name, *_ = test_results.keys()

        test_data = test_results[test_name]

        fragment = test_data[gene]['Amplicon']

        return fragment

    sequence = load_fragment(strain, gene, jsondir)

    # Leave at least a 50-mer on each end
    pivot = random.randint(50, len(sequence) - 50)
    side = random.choice((0, 1))

    halves = sequence[:pivot], sequence[:pivot]

    return halves[side]

def create_dummy_jsons(truncations: TRUNCATIONS, tempdir: Path) -> None:

    for strain in truncations:

        genes = {gene: {'Amplicon': truncations[strain][gene]}
                 for gene in truncations[strain]}

        data = {'Results': [{'TestResults': {'dummy': genes}}]}

        json_path = (tempdir / strain).with_suffix('.json')

        with json_path.open('w') as out:
            json.dump(data, out)


def recover_simulated(strain: str, gene: str, calls: pd.DataFrame,
                      distances: np.matrix, genes: Path, jsondir: Path,
                      replicates: int, seed: int) -> Union[int, str]:

    probs = crowbar.recover_allele(strain, gene, calls, distances, genes,
                                   jsondir, replicates, seed)


    most_likely_allele = max(probs, key=lambda x: probs[x])

    return most_likely_allele


def simulate_recovery(truncation_probability: float,
                      missing_probability: float, calls: pd.DataFrame,
                      jsondir: Path, genes: Path, tempdir: Path,
                      seed: int, replicates: int, cores: int):


    error_calls, truncations = random_errors(truncation_probability,
                                             missing_probability,
                                             calls, jsondir, seed)

    create_dummy_jsons(truncations, tempdir)

    distances = crowbar.hamming_distance_matrix(None, error_calls, cores)

    recover = partial(recover_simulated, calls=error_calls,
                      distances=distances, genes=genes, jsondir=tempdir,
                      seed=seed, replicates=replicates)

    result_futures = {}

    with ProcessPoolExecutor(max_workers=cores) as ppe:

        for strain in error_calls.index:

            for gene in error_calls.columns:

                if error_calls.loc[strain, gene] > 0:
                    continue

                recovered = ppe.submit(recover, strain, gene)

                result = {'recovered': recovered,
                          'original': calls.loc[strain, gene]}

                try:
                    result_futures[strain][gene] = result

                except KeyError:
                    result_futures[strain] = {}
                    result_futures[strain][gene] = result

    return result_futures

def main():

    args = arguments()

    assert args.tempdir != args.jsondir, 'tempdir cannot equal jsondir'

    recover(args.trunc_prob, args.miss_prob, args.distances, args.calls,
            args.jsons, args.genes, args.tempdir, args.seed, args.replicates,
            args.cores)

if __name__ == '__main__':
    main()
