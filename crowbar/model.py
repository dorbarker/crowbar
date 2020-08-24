import argparse
import shutil
import collections
import functools
import itertools
import json
import logging
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Dict, Tuple
import pandas as pd

from . import wallace
Abundance = Dict[str, Dict[str, int]]
Triplets = Dict[str, Dict[str, Tuple[str, float]]]


def build_model(calls_path: Path, alleles_dir: Path,
                model_path: Path, cores: int):

    logging.info('Creating model at `%s`', model_path)

    calls = load_calls(calls_path)

    model_path.mkdir(exist_ok=True, parents=True)

    save_calls(calls, model_path)

    save_known_alleles(alleles_dir, model_path)

    triplets = generate_copartitioning_triplets(calls, cores)

    abundances = calculate_abundances(calls)

    write_json(triplets, 'triplets.json', model_path)
    write_json(abundances, 'abundances.json', model_path)


def load_calls(calls_path: Path) -> pd.DataFrame:

    logging.info('Loading allele calls from `%s`', calls_path)

    raw_calls = pd.read_csv(calls_path, sep=',', index_col=0)

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

    namesA, namesB, valuesA, valuesB = make_mem_efficient_gene_arrays(calls)

    chunksize = math.ceil((len(calls.columns) ** 2) / (cores))

    with ProcessPoolExecutor(max_workers=cores) as ppe:

        results = ppe.map(wallace.bidirectional_adj_wallace, valuesA, valuesB,
                          chunksize=chunksize)

    for geneA, geneB, adj_wallace_value in zip(namesA, namesB, results):

        try:
            wallaces[geneA][geneB] = adj_wallace_value.a_b

        except KeyError:
            wallaces[geneA] = {geneB: adj_wallace_value.a_b}

        try:
            wallaces[geneB][geneA] = adj_wallace_value.b_a

        except KeyError:
            wallaces[geneB] = {geneA: adj_wallace_value.b_a}

    triplets = {}

    for geneA in wallaces:

        values = [(geneB, wallaces[geneB][geneA]) for geneB in wallaces[geneA]]

        best, second, *_ = sorted(values, key=lambda x: x[1], reverse=True)

        triplets[geneA] = {'best': best, 'second': second}

    return triplets


def make_mem_efficient_gene_arrays(calls: pd.DataFrame):

    names, values = zip(*calls.items())

    namesA, namesB = zip(*itertools.combinations(names, r=2))

    valuesA, valuesB = zip(*itertools.combinations(values, r=2))

    return namesA, namesB, valuesA, valuesB


def save_calls(calls: pd.DataFrame, model_path: Path) -> None:

    calls_path = model_path / 'calls.csv'

    logging.info('Save calls to %s', calls_path)
    calls.to_csv(calls_path)


def save_known_alleles(alleles_dir: Path, model_path: Path) -> None:

    if not alleles_dir.is_dir():
        logging.error('%s is not a directory', alleles_dir)
        sys.exit(1)

    model_alleles_dir = model_path / 'alleles'
    model_alleles_dir.mkdir(parents=True)

    # will copy all files in alleles_dir, as FASTA extensions aren't standard
    for fasta in alleles_dir.glob('*'):
        if fasta.is_file():
            dst = model_alleles_dir / fasta.name
            logging.info('Copying `%s` to `%s`', fasta, dst)
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

