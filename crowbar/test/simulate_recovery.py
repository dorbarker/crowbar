import sys
import os
import argparse
import json
import random
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from functools import partial
from typing import Dict, Tuple, Union
import pandas as pd

up = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
sys.path.append(up)

import crowbar  # main script

# Complex types
Truncations = Dict[str, str]
SimulationResults = Dict[str, Dict[str, Union[str, float]]]
PathTable = Dict[str, Path]

def arguments():


    parser = argparse.ArgumentParser()

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

    parser.add_argument('--test-jsons',
                        type=Path,
                        required=True,
                        help='Directory containing FSAC-format JSONs')

    parser.add_argument('--outdir',
                        type=Path,
                        required=True,
                        help='Output path')

    parser.add_argument('--model',
                        type=Path,
                        required=True,
                        help='Path to pre-trained model')

    parser.add_argument('-j', '--cores',
                        type=int,
                        default=1,
                        help='Number of CPU cores to use [1]')

    args = parser.parse_args()

    return args


def modify_row(strain_profile: pd.Series, trunc_count: int, miss_count: int,
               paths: PathTable) -> Tuple[pd.Series, Truncations]:

    to_truncate = random.sample(strain_profile.index, k=trunc_count)
    to_vanish = random.sample(strain_profile.index, k=miss_count)

    strain_profile[to_truncate] = -1
    strain_profile[to_vanish] = 0

    truncs = {gene: truncate(strain_profile.name, gene, paths)
              for gene in strain_profile[to_truncate].index}

    return strain_profile, truncs


def truncate(strain: str, gene: str, paths: PathTable) -> str:
    """Loads a partial nucleotide alignment from fsac-generated JSON
    output.

    Currently, this function assumes that there is only one cgMLST test per
    JSON file.

    :param strain:  The basename of the query strain
    :param gene:    The current query gene
    :param paths:   Dictionary of relevant paths
    :return:        A partial DNA sequence of at least length 50
    """

    jsonpath = paths['test'].joinpath(strain).with_suffix('.json')

    with jsonpath.open('r') as json_obj:
        data = json.load(json_obj)

    sequence = data[gene]['SubjAln']

    # Leave at least a 50-mer on each end
    pivot = random.randint(50, len(sequence) - 50)
    side = random.choice((0, 1))

    halves = sequence[:pivot], sequence[pivot:]

    return halves[side]


def create_dummy_jsons(strain: str, truncations: Truncations,
                       paths: PathTable) -> None:
    """Creates simplified JSONs containing sequence data for artificially
    truncated genes. These JSONs are compatible with FSAC JSONs.

    Because the recovery algorithm only refers to the JSONs in cases where
    partial sequence data are available, only genes with artificial truncations
    are written.

    :param strain:      Basename of the genome
    :param truncations: A Dict mapping partial sequence data to gene names
    :param paths:       Dictionary of relevant paths
    :return:            Path to the simulated FSAC JSON
    """

    genes = {gene: {'SubjAln': seq} for gene, seq in truncations.items()}

    sim_json_path = paths['simulated'].joinpath(strain).with_suffix('.json')

    with sim_json_path.open('w') as out:
        json.dump(genes, out)


def _simulate_recovery(genome_path: Path, trunc_count: int, miss_count: int,
                       paths: PathTable) -> Tuple[str, pd.Series]:
    """Simulates recovery of missing or truncated alleles by synthetically
    introducing these errors into error-free allele calls.

    After generating these synthetic errors and saves them into FSAC-formatted
    JSON files, it uses crowbar to attempt to recover these allele calls.

    :param trunc_count: Number of truncations to insert in each profile
    :param miss_count:  Number of missing alleles to insert in each profile
    :param paths:       Dictionary of relevant paths
                        JSONs created by this function
    :return:            Dict relating each strain name to its modified profile
    """


    strain_name = genome_path.stem

    strain_profile = crowbar.load_genome(genome_path)

    modified_profile, truncations = modify_row(strain_profile, trunc_count,
                                               miss_count, paths)

    create_dummy_jsons(strain_name, truncations, paths)

    repaired_calls, probabilities = crowbar.recover(modified_profile,
                                                    paths['simulated'],
                                                    paths['model'])

    crowbar.write_results(repaired_calls, probabilities, paths['recovered'])

    return strain_name, modified_profile


def simulate_recovery(trunc_count: int, miss_count: int, paths: PathTable,
                      cores: int) -> Dict[str, pd.Series]:


    jsons = paths['test'].glob('*.json')

    sim_recov = partial(_simulate_recovery,
                        trunc_count=trunc_count,
                        miss_count=miss_count,
                        paths=paths)

    with ProcessPoolExecutor(cores) as ppe:

        profiles = ppe.map(sim_recov, jsons)

    strain_profiles = dict(profiles)

    return strain_profiles


def compare_to_known(strain_profiles: Dict[str, pd.Series],
                     paths: PathTable) -> SimulationResults:
    # Compare maximum probabiltiy from JSON to known result
    # {gene: {allele: probability}}

    results = {}

    simulated_jsons = paths['recovered'].glob('*.json')

    for simulated_json in simulated_jsons:

        known_json = paths['test'] / simulated_json.name

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
                             'error_type': strain_profiles[strain][gene],
                             'correct': likeliest_allele == actual_allele,
                             'most_likely': likeliest_allele,
                             'most_likely_prob': alleles[likeliest_allele],
                             'second_allele': second_allele,
                             'second_prob': alleles[second_allele]}

    return results


def summarize_results(results: SimulationResults, paths: PathTable):

    df_results = pd.DataFrame(results).T

    df_results.to_csv(paths['report'], sep='\t', index=False)


def create_path_table(parent: Path, test_jsons: Path, model: Path) -> PathTable:
    """Helper function for ensuring that output paths are handled in a
    consistent manner.

    :param _parent:     The parent directory for all program output
    :param test_jsons:  Directory containing the query JSONs
    :param model:       Path to the pre-trained recovery model
    :return:            Dict containing output paths
    """

    paths = {
        'test':         test_jsons,
        'model':        model,
        'simulated':    parent / 'simulated_truncations',
        'recovered':    parent / 'recovered',
        'report':       parent / 'report.txt'
    }

    return paths


def main():

    args = arguments()

    paths = create_path_table(args.outdir, args.test_jsons, args.model)

    modified_profiles = simulate_recovery(args.trunc_count, args.miss_count,
                                          paths, args.cores)

    results = compare_to_known(modified_profiles, paths)

    summarize_results(results, paths)


if __name__ == '__main__':
    main()
