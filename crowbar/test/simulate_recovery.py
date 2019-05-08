# regular imports
import sys
import os
import argparse
import json
import random
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path
from typing import Dict, Optional, Tuple, Union
# 3rd Party imports
import pandas as pd
import numpy as np

up = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
sys.path.append(up)

import crowbar  # main script
from shared import logtime, hamming_distance_matrix

# Complex types
TRUNCATIONS = Dict[str, Dict[str, str]]
SimulationResults = Dict[str, Dict[str, Union[str, float]]]

def arguments():

    def get_temp_dir():

        return Path(os.environ.get('XDG_RUNTIME_DIR') or '/tmp')

    parser = argparse.ArgumentParser()

    parser.add_argument('--cores',
                        type=int,
                        default=1,
                        help='Number of CPU cores to use [1]')

    parser.add_argument('--truncation-probability',
                        type=float,
                        default=0.0,
                        dest='trunc_prob',
                        help='Uniform probability that any given \
                              locus will be truncated [0.0]')

    parser.add_argument('--missing-probability',
                        type=float,
                        default=0.0,
                        dest='miss_prob',
                        help='Uniform probability that any given \
                              locus will be rendered missing [0.0]')

    parser.add_argument('--tempdir',
                        type=Path,
                        default=get_temp_dir(),
                        help='Directory for ephemeral working files')

    parser.add_argument('--output',
                        type=Path,
                        required=True,
                        help='Output path')

    parser.add_argument('--model',
                        type=Path,
                        required=True,
                        help='Path to pre-trained model')

    args = parser.parse_args()

    model_jsons = args.model / 'jsons'

    if args.tempdir == model_jsons:

        msg = '--tempdir [{}] cannot be the JSON directory of the model [{}]'
        sys.exit(msg.format(args.tempdir, model_jsons))

    return args


def modify_row(strain_profile: pd.Series, trunc_count: int, miss_count: int,
               jsondir: Path) -> Tuple[pd.Series, Dict[str, str]]:

    to_truncate = random.sample(strain_profile.index, k=trunc_count)
    to_vanish = random.sample(strain_profile.index, k=miss_count)

    strain_profile[to_truncate] = -1
    strain_profile[to_vanish] = 0

    truncs = {gene: truncate(strain, gene, jsondir)
              for gene in row[to_truncate].index}

    return strain_profile, truncs


def truncate(strain: str, gene: str, jsondir: Path) -> str:
    """Loads a partial nucleotide alignment from fsac-generated JSON
    output.

    Currently, this function assumes that there is only one cgMLST test per
    JSON file.

    :param strain:  The basename of the query strain
    :param gene:    The current query gene
    :param jsondir: Path to the directory containing FSAC input JSONs

    :return:        A partial DNA sequence of at least length 50
    """

    jsonpath = (jsondir / strain).with_suffix('.json')

    with jsonpath.open('r') as json_obj:
        data = json.load(json_obj)

    sequence = data[gene]['SubjAln']

    # Leave at least a 50-mer on each end
    pivot = random.randint(50, len(sequence) - 50)
    side = random.choice((0, 1))

    halves = sequence[:pivot], sequence[pivot:]

    return halves[side]


@logtime('Creating dummy JSONs')
def create_dummy_jsons(strain: str, truncations: TRUNCATIONS,
                       tempdir: Path) -> Path:

    genes = {gene: {'SubjAln': seq} for gene, seq in truncations.items()}

    temp_json_path = (tempdir / strain).with_suffix('.json')

    with temp_json_path.open('w') as out:
        json.dump(genes, out)

    return temp_json_path


def simulate_recovery(trunc_count: int, miss_count: int, json_dir: Path,
                      temp_dir: Path, outdir: Path, model_path: Path):
    """Simulates recovery of missing or truncated alleles by synthetically
    introducing these errors into error-free allele calls.

    After generating these synthetic errors and saves them into FSAC-formatted
    JSON files, it uses crowbar to attempt to recover these allele calls.

    :param trunc_count: Number of truncations to insert in each profile
    :param miss_count:  Number of missing alleles to insert in each profile
    :param json_dir:    Location of the directory containing FSAC JSONs
    :param temp_dir:    Location of the directory containing temporary dummy
                        JSONs created by this function
    :param outdir:      Location to write crowbar recovery JSONs
    :param model_path:  Location of a pre-trained crowbar model
    """
    jsons = json_dir.glob('*.json')

    for genome_path in jsons:

        strain_name = genome_path.stem

        strain_profile = crowbar.load_genome(genome_path)

        modified_profile, truncation = modify_row(strain_profile, trunc_count,
                                                  miss_count, json_dir)

        temp_json = create_dummy_jsons(strain_name, truncations, tempdir)

        repaired_calls, probabilities = crowbar.recover(modified_profile,
                                                        temp_json, model_path)

        crowbar.write_results(repaired_calls, probabilities, outdir)


def compare_to_known(simulation_outdir: Path,
                     known_jsondir: Path,
                     strain_profile: pd.Series) -> SimulationResults:
    # Compare maximum probabiltiy from JSON to known result
    # {gene: {allele: probability}}

    results = {}

    simulated_jsons = simulation_outdir.glob('*.json')

    for simulated_json in simulated_jsons:

        known_json = known_jsondir / simulated_json.name

        strain = simulated_json.stem

        with simulated_json.open('r') as s, known_json.open('r') as k:

            simulated = json.load(s)

            known = json.load(k)

        for gene in simulated:

            alleles = simulated[gene]

            likeliest_allele, second_allele = sorted(alleles,
                                                     key=lambda x: alleles[x])

            actual_allele = known[gene]['MarkerMatch']

            results[gene] = {'actual': actual_allele,
                             'error_type': strain_profile[gene],
                             'correct': likeliest_allele == actual_allele,
                             'most_likely': likeliest_allele,
                             'most_likely_prob': alleles[likeliest_allele],
                             'second_allele': second_allele,
                             'second_prob': alleles[second_allele]}

    return results


def summarize_results(results: SimulationResults, result_out: Path):

    df_results = pd.DataFrame(results).T

    df_results.to_csv(result_out, sep='\t', index=False)


def main():

    args = arguments()


    summarize_results(results, args.output)

if __name__ == '__main__':
    main()
