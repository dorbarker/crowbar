import argparse
import shutil
import collections
import functools
import itertools
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Dict, Tuple
import pandas as pd

import wallace
Abundance = Dict[str, Dict[str, int]]
Triplets = Dict[str, Dict[str, Tuple[str, float]]]

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('--calls',
                        type=Path,
                        required=True,
                        help='Path to allele calls table')

    parser.add_argument('--alleles-dir',
                        type=Path,
                        required=True,
                        help='Path to directory containing gene FASTAs')

    parser.add_argument('--output',
                        type=Path,
                        required=True,
                        help='Directory containing model')

    parser.add_argument('--cores',
                        type=int,
                        default=1,
                        help='Number of CPU cores to use [1]')

    return parser.parse_args()


def main():

    args = arguments()

    calls = load_calls(args.calls)

    build_model(calls, args.alleles_dir, args.output, args.cores)


def build_model(calls: pd.DataFrame, alleles_dir: Path,
                model_path: Path, cores: int):

    model_path.mkdir(exist_ok=True, parents=True)

    save_calls(calls, model_path)

    save_known_alleles(alleles_dir, model_path)

    triplets = generate_copartitioning_triplets(calls, cores)

    abundances = calculate_abundances(calls)

    write_json(triplets, 'triplets.json', model_path)
    write_json(abundances, 'abundances.json', model_path)


def load_calls(calls_path: Path) -> pd.DataFrame:

    raw_calls = pd.read_csv(calls_path, index_col=0, dtype=int)

    calls = raw_calls.loc[[not any(v < 1) for i, v in raw_calls.iterrows()]]

    return calls


def generate_copartitioning_triplets(calls: pd.DataFrame, cores: int) -> Triplets:
    """For each gene, find the other two genes that partition the population
    most similarly to that gene. Uses the Adjusted Wallace Coefficient to
    determine the most similar genes.

    :param calls:   DataFrame containing allele calls
    :return:        Dictionary of type Triplets
    """

    wallaces = {}

    geneAs, geneBs = zip(*itertools.permutations(calls.columns, r=2))

    adj_wallace = functools.partial(wallace.adj_wallace, calls=calls)

    with ProcessPoolExecutor(max_workers=cores) as ppe:

        results = ppe.map(adj_wallace, geneAs, geneBs, chunksize=100)

    for geneA, geneB, adj_wallace_value in zip(geneAs, geneBs, results):

        try:
            wallaces[geneA][geneB] = adj_wallace_value

        except KeyError:
            wallaces[geneA] = {geneB: adj_wallace_value}

    triplets = {}

    for geneA in wallaces:

        values = [(geneB, wallaces[geneB][geneA]) for geneB in wallaces[geneA]]

        best, second, *_ = sorted(values, key=lambda x: x[1], reverse=True)

        triplets[geneA] = {'best': best, 'second': second}

    return triplets


def save_calls(calls: pd.DataFrame, model_path: Path) -> None:

    calls_path = model_path / 'calls.csv'

    calls.to_csv(calls_path)


def save_known_alleles(alleles_dir: Path, model_path: Path) -> None:

    if not alleles_dir.is_dir():
        msg = f"{alleles_dir} is not a directory."
        raise NotADirectoryError(msg)

    model_alleles_dir = model_path / 'alleles'
    model_alleles_dir.mkdir(parents=True)

    # will copy all files in alleles_dir, as FASTA extensions aren't standard
    for fasta in alleles_dir.glob('*'):
        if fasta.is_file():
            dst = model_alleles_dir / fasta.name
            shutil.copy(fasta, dst)


def calculate_abundances(calls: pd.DataFrame) -> Abundance:
    """Calculate the absolute abundances of non-missing alleles for each gene
    in the calls table. Also treat add an unknown allele represented by '?'
    with a frequency of 1.

    :param calls:   A dataframe of allele calls with each
                    row a strain and each column a gene.

    :return:        A dictionary of allele frequencies for each gene.
    """

    abundances = {}

    for gene, alleles in calls.iteritems():

        present = alleles[alleles > 0]

        inc_unknown = [str(x) for x in present] + ['?']

        frequencies = collections.Counter(inc_unknown)

        abundances[gene] = dict(frequencies)  # simplify for later JSON

    return abundances


def write_json(values: Dict, name: str, model_path: Path) -> None:
    """Saves dictionary data as a JSON-formatted file."""

    output_path = model_path / name

    with output_path.open('w') as f:
        json.dump(values, f, indent=4)


def calculate_distance_matrix(calls_path: Path, model_path: Path) -> None:
    """Calculates a hamming distance matrix using the external `hamming` tool.
    https://gitlab.com/dorbarker/hamming

    cargo install --git https://gitlab.com/dorbarker/hamming.git

    :param calls_path: Path to the input allele calls
    :param model_path: Directory containing the model
    """

    cmd = ('hamming', '--input', str(calls_path), '--output', str(model_path))

    subprocess.run(cmd, check=True)


if __name__ == '__main__':
    main()
