import argparse
import shutil
import collections
import itertools
import json
from pathlib import Path
from typing import Dict
import pandas as pd

import wallace
Abundance = Dict[str, Dict[str, int]]


def arguments():

    parser = argparse.ArgumentParser()

    return parser.parse_args()


def main():
    ...


def build_model():
    ...


def load_calls(calls_path: Path) -> pd.DataFrame:

    raw_calls = pd.read_csv(calls_path, index_col=0)

    calls = raw_calls.loc[[not any(v < 1) for i, v in raw_calls.iterrows()]]

    return calls


def reorder_calls():

    wallaces = {}

    pairwise_calls = itertools.permutations(calls.iteritems(), r=2)

    for (geneA, allelesA), (geneB, allelesB) in pairwise_calls:

        adj_wallace_value = wallace.adj_wallace(allelesA, allelesB)

        try:
            wallaces[geneA][geneB] = adj_wallace_value

        except KeyError:
            wallaces[geneA] = {geneB: adj_wallace_value}






def save_calls(calls: pd.DataFrame, model_path: Path) -> None:

    calls.to_csv(model_path)


def save_known_alleles(alleles_dir: Path, model_path: Path) -> None:

    if not alleles_dir.is_dir():
        msg = f"{alleles_dir} is not a directory."
        raise NotADirectoryError(msg)

    model_alleles_dir = model_path / 'alleles'
    model_alleles_dir.mkdir(parents=True)

    # will copy all files in alleles_dir, as FASTA extensions aren't standard
    for fasta in alleles_dir.glob('*'):
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


def save_abundances(abundances: Abundance, model_path: Path) -> None:
    """Saves allele abundance data as a JSON-formatted file."""

    abundance_path = model_path / 'allele_abundances.json'

    with abundance_path.open('w') as f:
        json.dump(abundances, f, indent=4)


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
