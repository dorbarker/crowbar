import argparse
import itertools
import json
import re
import math
import operator
import subprocess
import sys
import warnings
from collections import Counter, namedtuple
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from io import StringIO
from pathlib import Path
from statistics import mean
from typing import Dict, List, Optional, Set, Tuple, Union

# Third-party imports
import pandas as pd
import numpy as np
from Bio import SeqIO

# Local imports
from shared import user_msg, logtime, hamming_distance_matrix

# Complex type constants
TRIPLET_COUNTS = Dict[int, Dict[int, Dict[int, int]]]
TRIPLET_PROBS = Dict[int, Dict[int, Dict[int, float]]]
NUMERIC = Union[int, float]
ABUNDANCE = Dict[Union[str, int], float]

Neighbour = namedtuple('Neighbour', ('indices', 'alleles', 'similarity'))

NeighbourAlleles = List[Tuple[str, float]]
AlleleProb = Dict[str, float]

def arguments():
    """Gather command line arguments for crowbar.py"""

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input',
                        type=Path,
                        required=True,
                        help='Path to input fsac-formatted JSON')

    parser.add_argument('-o', '--output',
                        type=Path,
                        required=True,
                        help='Path to output directory')

    parser.add_argument('--model',
                        type=Path,
                        required=True,
                        help='Path to recovery model')

    return parser.parse_args()


def nearest_neighbour(strain_profile: pd.Series, gene: str,
                      model_calls: pd.DataFrame) -> NeighbourAlleles:
    """Finds the nearest neighbour(s) to `strain`, and returns a Neighbour
    namedtuple containing the row indices in `calls` of its closest relatives,
    the allele(s) present at `gene`, and their percent Hamming similarity.
    """

    def percent_shared(strain1: pd.Series, strain2: pd.Series) -> float:
        """Returns the percent similarity of two strains based on the Hamming
        distance of non-missing loci.
        """

        shared = [a and b for a, b in zip(strain1 > 0, strain2 > 0)]

        return sum(strain1[shared] == strain2[shared]) / sum(shared)


    similarities = [(row_strain, percent_shared(strain_profile, row_profile))
                    for row_strain, row_profile in model_calls.iterrows()]

    sorted_similarities = reversed(sorted(similarities, key=lambda x: x[1]))

    max_similarity = sorted_similarities[0][1]

    neighbours = itertools.takewhile(lambda x: x[1] == max_similarity,
                                     sorted_similarities)

    neighbouring_alleles = [(model_calls.loc[row_strain, gene], similarity)
                            for row_strain, similarity in neighbours]

    return neighbouring_alleles


def neighbour_allele_probabilities(neighbouring_alleles: NeighbourAlleles,
                                   abundances: AlleleProb) -> AlleleProb:
    ...

    neighbours, similarities = zip(*neighbouring_alleles)
    similarity, *_ = similarities

    allele_counts = Counter(neighbours)

    for allele, abundance in abundances.items():

        if allele in neighbours:

            count = allele_counts[allele]

            combined_probs[allele] = similarity * (1 - (abundance ** count))

        else:

            combined_probs[allele] = (1 - similarity) * abundance

    return combined_probs


def closest_relative_alleles(strain_profile: pd.Series, gene: str,
                             model_path: Path,
                             abundances: AlleleProb) -> AlleleProb:

    calls = pd.read_csv(model_path / 'calls.csv')

    neighbour_alleles = nearest_neighbour(strain_profile, gene, calls)

    return neighbour_allele_probabilities(neighbour_alleles, abundances)

def flank_linkage(strain_profile: pd.Series, gene: str, model_path: Path,
                  abundances: Dict) -> AlleleProb:

    calls_path = model_path / calls.csv
    calls = pd.read_csv(calls_path)

    triplet_path = model_path / 'triplets.json'
    with triplet_path.open('r') as f:
        data = json.load(f)

    best_predictor_gene = data[gene]['best']
    second_predictor_gene = data[gene]['second']

    best_calls = calls[best_predictor_gene]
    second_calls = calls[second_predictor_gene]
    query_calls = calls[gene]

    possible = np.array(abundances.keys())

    strain_best_predictor = strain_profile[best_predictor_gene]
    strain_second_predictor = strain_profile[second_predictor_gene]

    predictions = (best_calls == strain_best_predictor)     & \
                  (second_calls == strain_second_predictor) & \
                  (np.isin(query_calls, possible))

    triplet_probabilities = {}

    for hypothesis in possible:

        if hypothesis not in query_calls:
            # If hypothesis is *never* observed with these predictors,
            # it's not impossible - merely unlikely.
            #
            # Use the probability of an unobserved allele
            # as the base probability.

            triplet_probabilities[hypothesis] = abundances['?']

        else:

            is_h = (calls[gene] == hypothesis)

            predictors_given_h = (has_predictors & is_h).sum() / is_h.sum()

            triplet_probabilities[hypothesis] = predictors_given_h

    return triplet_probabilities


def partial_sequence_match(gene: str, model_path: Path, jsonpath: Path) -> Set[str]:
    """Attempts to use partial sequence data to exclude possible alleles."""

    genes = model_path / 'alleles'

    with jsonpath.open('r') as json_obj:
        data = json.load(json_obj)

    fragment = data[gene]['SubjAln']

    fragment_pattern = re.compile('(^{seq})|({seq}$)'.format(seq=fragment))

    glob_pattern = '*{}.f*'.format(gene)
    gene_file, *_ = genes.glob(glob_pattern)

    with gene_file.open('r') as fasta:

        matches = set(rec.id for rec in SeqIO.parse(fasta, 'fasta')
                      if re.search(fragment_pattern, str(rec.seq)))

    return matches


def allele_abundances(gene: str, partial_matches: Set[str],
                      model_path: Path) -> AlleleProb:
    """Calculates the abundances of alleles for a given gene and attempts to
    determine the probability that the next observation will be a new allele.
    """

    abundance_path = model_path / 'abundances.json'

    with abundance_path.open('r') as f:
        data = json.load(f)

    model_abundances = data[gene]

    possible_abundance = {allele: count
                          for allele, count in model_abundances.items()
                          if allele in partial_matches or allele == '?'}

    total_observations = sum(possible_abundance.values())

    abundance = {allele: count / total_observations
                 for allele, count in possible_abundance}

    return abundance


def bayes(abundances: AlleleProb, triplets: AlleleProb,
          neighbours: AlleleProb) -> AlleleProb:

    def bayes_theorem(h: Union[str, int]) -> float:
        """Implementation of Bayes' Theorem in which the hypothesis being
        tested (h) is a given allele, and the lines of evidence are
        linkage disequilibrium between genes, and the similarity of the strain
        in question to its nearest neighbour.
        """

        neighbour_probability = neighbours[h]

        triplet_probability = triplets[h]

        e_h = triplet_probability * neighbour_probability

        return e_h

    likelihoods = {h: bayes_theorem(h) for h in adj_abundances}

    e = sum(likelihoods.values())

    return {h: ((likelihoods[h] * abundances[h]) / e) for h in abundances}


def load_genome(genome_calls_path: Path):

    # Load from fsac JSON
    with genome_calls_path.open('r') as f:
        data = json.load(f)

    genome_calls = {}
    for gene in data:

        if data[gene]['BlastResult'] is False:
            genome_calls[gene] = '0'

        elif data[gene]['IsContigTrucation']:
            genome_calls[gene] = '-1'

        else:
            genome_calls[gene] = data[gene]['MarkerMatch']

    return pd.Series(genome_calls, name=genome_calls_path.stem)


def recover(strain_profile: pd.Series, json_path,
            model_path: Path) -> Tuple[pd.Series, Dict[str, AlleleProb]]:


    repaired_calls = strain_profile.copy()

    all_gene_probabilities = {}

    missing = strain_profile < 1

    missing_genes = strain_profile[missing].index

    for gene in missing_genes:

        partial_matches = partial_sequence_match(gene, model_path, json_path)

        abundances = allele_abundances(gene, partial_matches, model_path)

        triplets = flank_linkage(strain_profile, gene, model_path, abundances)

        neighbours = closest_relative_alleles(strain_profile, gene,
                                              model_path, abundances)

        probabilities = bayes(abundances, triplets, neighbours)

        most_probable = max(probabilities, key=lambda x: probabilities[x])

        repaired_calls[gene] = most_probable

        all_gene_probabilities[gene] = probabilities

    return repaired_calls, all_gene_probabilities


def write_results(repaired_calls: pd.Series,
                  all_gene_probabilities: Dict[str, AlleleProb],
                  outdir: Path) -> None:

    output_name = outdir / repaired_calls.name

    pd.DataFrame(repaired_calls).T.to_csv(output_name.with_suffix('.csv'))

    with output_name.with_suffix('.json').open('w') as f:
        json.dump(all_gene_probabilities)


def main():
    """Main function. Gathers arguments and passes them to recover()"""

    args = arguments()

    strain_profile = load_genome(args.input)

    repaired, probabilities =  recover(strain_profile, args.input, args.model)

    write_results(repaired, probabilities, args.output)


if __name__ == '__main__':
    main()
