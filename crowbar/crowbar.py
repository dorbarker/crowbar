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

    parser.add_argument('--cores',
                        type=int,
                        default=1,
                        help='Number of CPU cores to use [1]')

    parser.add_argument('--replicates',
                        type=int,
                        default=100,
                        help='Replicates for the Monte Carlo estimation of \
                              the probability of finding a novel allele')

    parser.add_argument('--seed',
                        type=int,
                        default=1,
                        help='Seed to initialize the Monte Carlo simulation')

    parser.add_argument('--reference',
                        type=Path,
                        required=True,
                        help='Path to reference genome')

    parser.add_argument('--distances',
                        type=Path,
                        required=False,
                        help='Path to pre-calculated distance matrix')

    parser.add_argument('--output',
                        type=Path,
                        required=False,
                        help='JSON-formatted output destination [stdout]')

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


def nearest_neighbour(gene: str, strain_profile: pd.Series,
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


def flank_linkage(strain_profile: pd.Series, gene: str, model_path: Path,
                  abundances: Dict, calls: pd.DataFrame) -> AlleleProb:

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


def partial_sequence_match(gene: str, genes: Path, jsonpath: Path) -> Set[str]:
    """Attempts to use partial sequence data to exclude possible alleles."""


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


def allele_abundances(gene: str, partial_matches: Set[str], model_path: Path):
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


def bayes(adj_abundances, neighbour_probs, flanks) -> Dict[int, float]:

    def bayes_theorem(h: Union[str, int]) -> float:
        """Implementation of Bayes' Theorem in which the hypothesis being
        tested (h) is a given allele, and the lines of evidence are
        linkage disequilibrium between genes, and the similarity of the strain
        in question to its nearest neighbour.
        """

        neighbour_probability = neighbour_probs[h]

        flanks_given_h = flanks[h]

        e_h = flanks_given_h * neighbour_probability

        return e_h

    likelihoods = {h: bayes_theorem(h) for h in adj_abundances}

    e = sum(likelihoods.values())

    return {h: ((likelihoods[h] * adj_abundances[h]) / e)
            for h in adj_abundances}


def bayes(abundances, triplets, neighbours):
    ...


def recover_allele(strain: str, gene: str, calls: pd.DataFrame,
                   distances: np.matrix, genes: Path, jsondir: Path,
                   gene_abundances):
    """Attempts to recover the allele of a missing locus based on
    lines of evidence from allele abundance, linkage disequilibrium, the
    close relatives.

    If the locus is truncated, fragment matching can be used to eliminate
    possible alleles. If the gene is wholly missing, then all alleles are in
    contention.
    """

    if calls.loc[strain, gene] == -1:  # is truncation

        fragment_matches = partial_sequence_match(strain, gene, genes, jsondir)

    else:  # is wholly missing

        fragment_matches = set(calls[gene]) - {0, -1}

    adj_abundances = redistribute_allele_probability(gene_abundances[gene],
                                                     fragment_matches)

    neighbour = nearest_neighbour(gene, strain,  # fragment_matches,
                                  distances, calls)

    neighbour_probs = neighbour_similarities(neighbour, adj_abundances)

    flanks = {h: flank_linkage(strain, gene, h, gene_abundances, calls)
              for h in adj_abundances}

    probabilities = bayes(adj_abundances, neighbour_probs, flanks)

    return probabilities


def recover(callspath: Path, reference: Path, genes: Path, jsondir: Path,
            distance_path: Optional[Path], replicates: int, seed: int,
            cores: int):
    """Master function for crowBAR.

    Walks through a pandas DataFrame of allele calls and attempts to recover
    any truncated (-1) or absent (0) loci.
    """

    results = {}  # Dict[str, Dict[str, Dict[int, float]

    calls = order_on_reference(reference, genes,
                               pd.read_csv(callspath, index_col=0), cores)

    distances = hamming_distance_matrix(distance_path, calls, cores)

    gene_abundances = gene_allele_abundances(calls, replicates, seed, cores)

    total_bad = sum((calls < 1).sum())

    counter = 0

    for strain in calls.index:

        for gene in calls.columns:

            if calls.loc[strain, gene] < 1:  # truncated or missing

                probs = recover_allele(strain, gene, calls, distances,
                                       genes, jsondir, gene_abundances)

                try:
                    results[strain][gene] = probs

                except KeyError:

                    results[strain] = {}
                    results[strain][gene] = probs

                counter += 1

                user_msg(counter, '/', total_bad)

    return results


def write_output(results, outpath: Path) -> None:

    reformatted_results = {}

    for strain in results:
        reformatted_results[strain] = {}

        for gene in results[strain]:
            reformatted_results[strain][gene] = {}

            for allele in results[strain][gene]:

                allele_ = str(allele)  # JSON keys are converted to str anyway

                prob = results[strain][gene][allele]
                probability = -1 if np.isnan(prob) else float(prob)

                reformatted_results[strain][gene][allele_] = probability

    if outpath is None:

        json.dump(reformatted_results, sys.stdout, indent=4)

    else:

        with outpath.open('w') as out:
            json.dump(reformatted_results, out, indent=4)


def main():
    """Main function. Gathers arguments and passes them to recover()"""

    args = arguments()

    results = recover(args.calls, args.reference, args.genes, args.jsons,
                      args.distances, args.replicates, args.seed, args.cores)

    write_output(results, args.output)


if __name__ == '__main__':
    main()
