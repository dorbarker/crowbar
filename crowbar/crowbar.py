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
                      model_calls: pd.DataFrame) -> Neighbour:
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


def find_locus_in_reference(gene: Path, reference: Path) -> Tuple[str, int]:

    def find_locus_location(seq: str) -> int:
        """Executes a BLASTn search for a core gene against a reference genome.

        Returns the minimum of the start and stop locations of the gene
        so that gene calls can be ordered relative to the reference genome.
        """

        blast = ('blastn', '-task', 'megablast', '-subject', str(reference),
                 '-outfmt', '10')

        string_result = subprocess.run(blast, universal_newlines=True,
                                       check=True, input=seq,
                                       stdout=subprocess.PIPE)

        table_result = pd.read_table(StringIO(string_result.stdout),
                                     sep=',', header=None)

        # best row should be first
        # sstart, ssend = 8, 9
        start, stop = table_result.iloc[0, [8, 9]]

        return min(start, stop)

    gene_name = gene.stem

    with gene.open('r') as fasta:

        # Use just the first record
        rec = next(SeqIO.parse(fasta, 'fasta'))

        query = str(rec.seq)

        loc = find_locus_location(query)

    return gene_name, loc


def flank_linkage(strain: str, gene: str, hypothesis: int, gene_abundances,
                  calls: pd.DataFrame) -> Tuple[float, float]:
    """For three genes, (Left, Centre, Right), count how many observations of
    Centre were associated with the flanking genes Left and Right.
    """

    possible = np.array([a for a in gene_abundances[gene].keys() if a != '?'])

    columns = tuple(calls.columns)
    gene_loc = columns.index(gene)

    left_col = columns[gene_loc - 1]

    # Wrap around if the centre gene is already the rightmost column
    try:
        right_col = columns[gene_loc + 1]
    except IndexError:
        right_col = columns[0]

    flank_left, flank_right = calls[[left_col, right_col]].loc[strain]

    left = np.array(calls[left_col])
    rght = np.array(calls[right_col])
    cntr = np.array(calls[gene])

    has_flank = (left == flank_left) & \
                (rght == flank_right) & \
                (np.isin(cntr, possible))

    with warnings.catch_warnings():

        warnings.simplefilter('ignore')

        try:
            is_h = (calls[gene] == hypothesis)

            flanks_given_h = (has_flank & is_h).sum() / is_h.sum()

            if math.isnan(flanks_given_h):
                raise TypeError

        except TypeError:

            flanks_given_h = 0

    # If hypothesis is *never* observed with these flanks, it's not
    # impossible - merely unlikely. Use the probability of an unobserved allele
    # as the base probability
    return flanks_given_h or gene_abundances[gene]['?']


def partial_sequence_match(strain: str, gene: str, genes: Path,
                           jsondir: Path) -> Set[int]:
    """Attempts to use partial sequence data to exclude possible alleles."""

    def load_fragment() -> str:
        """Loads a partial nucleotide alignment from FSAC-generated JSON output.

        :return: A partial gene alignment as a str
        """

        jsonpath = (jsondir / strain).with_suffix('.json')

        with jsonpath.open('r') as json_obj:
            data = json.load(json_obj)

        return data[gene]['SubjAln']

    def fragment_match(seq: str) -> Set[int]:
        """Attempts to match partial sequence data to a known allele from
        a multifasta file. Matches are only attemped at the beginning and end
        of the gene.
        """

        fragment_pattern = re.compile('(^{seq})|({seq}$)'.format(seq=seq))

        glob_pattern = '*{}.f*'.format(gene)
        gene_file, *_ = genes.glob(glob_pattern)

        with gene_file.open('r') as fasta:

            result = set(int(rec.id) for rec in SeqIO.parse(fasta, 'fasta')
                         if re.search(fragment_pattern, str(rec.seq)))

        return result

    fragment = load_fragment()

    matches = fragment_match(fragment)

    return matches


def allele_abundances(gene: str, calls: pd.DataFrame, replicates: int = 1000,
                      seed: int = 1) -> Tuple[str, ABUNDANCE]:
    """Calculates the abundances of alleles for a given gene and attempts to
    determine the probability that the next observation will be a new allele.
    """

    def account_for_unknown(freq: int, observations: int,
                            unknown: float) -> float:

        return (freq / observations) - unknown

    observed_alleles = richness_estimate.Population(calls[gene])

    n_alleles = len(observed_alleles.abundance)

    discoveries = observed_alleles.monte_carlo(replicates, seed)

    last_percentile = max(1, int(0.01 * len(discoveries)))

    last_percentile_discovery_rate = mean(discoveries[-last_percentile:])

    to_subtract_from_known = last_percentile_discovery_rate / n_alleles

    account_for_unknown_ = partial(account_for_unknown,
                                   observations=len(observed_alleles),
                                   unknown=to_subtract_from_known)

    abundance = {allele: account_for_unknown_(frequency)
                 for allele, frequency
                 in observed_alleles.abundance.items()}

    abundance['?'] = last_percentile_discovery_rate

    return gene, abundance


def redistribute_allele_probability(abundances, fragment_matches):
    """Filters alleles ruled out by fragment matching and recalculates
    abundance proportions.
    """

    filtered_abundances = {k: v
                           for k, v in abundances.items()
                           if k in fragment_matches or k == '?'}

    total_remaining_probability = sum(filtered_abundances.values())

    adjusted_abundances = {k: v / total_remaining_probability
                           for k, v in filtered_abundances.items()}

    return adjusted_abundances


def neighbour_similarities(neighbour: Neighbour,
                           abundances: Dict[Union[int, str], float]):
    """Use weighting for multiple observations of the same neightbour allele"""

    allele_proportions = Counter(neighbour.alleles)

    combined_probs = {}

    for k, abund in abundances.items():

        if k in allele_proportions:

            count = allele_proportions[k]

            combined_probs[k] = neighbour.similarity * (1 - (abund ** count))

        else:

            combined_probs[k] = (1 - neighbour.similarity) * abund

    return combined_probs


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
