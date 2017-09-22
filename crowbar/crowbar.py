import argparse
import json
import re
import math
import operator
import subprocess
import warnings
from collections import Counter, namedtuple
from concurrent.futures import ProcessPoolExecutor
from functools import partial, reduce
from io import StringIO
from pathlib import Path
from statistics import mean
from typing import Callable, Dict, List, Optional, Sequence, Set, Tuple, Union

# Third-party imports
import pandas as pd
import numpy as np
from Bio import SeqIO

# Local imports
import richness_estimate
from shared import user_msg, logtime

# Complex type constants
TRIPLET_COUNTS = Dict[int, Dict[int, Dict[int, int]]]
TRIPLET_PROBS = Dict[int, Dict[int, Dict[int, float]]]
TRIPLET_COUNT_CURRENT = Tuple[TRIPLET_COUNTS, Dict[int, int]]
NUMERIC = Union[int, float]

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

@logtime('Getting allele abundances')
def gene_allele_abundances(calls: pd.DataFrame, replicates: int,
                           seed: int, cores: int) -> Dict[str, float]:
    """Get the allele probabilities, including the those of a novel allele,
    for every gene in `calls`.

    Runs in parallel.
    """

    with ProcessPoolExecutor(max_workers=cores) as ppe:

        get_allele_abundances = partial(allele_abundances,
                                        calls=calls,
                                        replicates=replicates,
                                        seed=seed)

        abundances = dict(ppe.map(get_allele_abundances, calls.columns))

    return abundances


def row_distance(idx: int, row: Tuple[str, pd.Series],
                 calls: pd.DataFrame) -> Dict[int, int]:
    """Returns the Hamming distance of non-missing alleles between two strains.

    Results are returned as a dictionary for distances between the query strain
    (strain1) and all the subject strain (each strain2).

    Called by dist_gene()
    """

    strain1 = row[1]

    def non_missing_hamming(j):
        """Returns the distance between two strains, considering only loci
        which are not missing in either individual.
        """

        strain2 = calls.iloc[j]

        return sum([a > 0 and b > 0 and a != b
                    for a, b in zip(strain1, strain2)])

    return {j: non_missing_hamming(j) for j in range(idx + 1, len(calls))}


@logtime('Calculating distance matrix')
def hamming_distance_matrix(distance_path: Optional[Path], calls: pd.DataFrame,
                            cores: int) -> np.matrix:

    """Returns a Hamming distance matrix of pairwise strain distances.

    First attempts to load a pre-calculated matrix from `distance_path and
    returns that if possible. If there is no distance matrix located at that
    path, it calculates one and saves it there. If no path is provided,
    calculate the distance matrix and return it without saving to disk.
    """

    def dist_gene(calls: pd.DataFrame, cores: int) -> np.matrix:
        """Returns a Hamming distance matrix of pairwise strain distances."""

        n_row = len(calls)

        dist_mat = np.matrix([np.zeros(n_row)
                              for _ in range(n_row)],
                              dtype=int)

        with ProcessPoolExecutor(max_workers=cores) as ppe:
            futures = {i: ppe.submit(row_distance, i, row, calls)
                       for i, row in enumerate(calls.iterrows())}

        results = {i: j.result() for i, j in futures.items()}

        for i, js in results.items():

            for j in js:
                dist_mat[i, j] = dist_mat[j, i] = results[i][j]

        return dist_mat


    try:

        distances = np.matrix(pd.read_csv(distance_path,
                                          header=0, index_col=0))

    except FileNotFoundError:

        distances = dist_gene(calls, cores)

        pd.DataFrame(distances,
                     index=calls.index,
                     columns=calls.index).to_csv(distance_path)

    except ValueError:

        distances = dist_gene(calls, cores)

    return distances

Neighbour = namedtuple('Neighbour', ('indices', 'alleles', 'similarity'))
def nearest_neighbour(gene: str, strain: str, included_fragments: Set[int],
                      distances: np.matrix, calls: pd.DataFrame) -> Neighbour:
    """Finds the nearest neighbour(s) to `strain`, and returns a Neighbour
    namedtuple containing the row indices in `calls` of its closest relatives,
    the allele(s) present at `gene`, and their percent Hamming similarity.
    """

    def percent_shared(strain1: pd.Series, strain2: pd.Series) -> float:
        """Returns the percent similarity of two strains based on the Hamming
        distance of non-missing loci.
        """

        shared = [a and b for a, b in zip(strain1 > 0, strain2 > 0)]

        return sum(strain1[shared] == strain2[shared]) / len(shared)


    def which(data: Sequence, compared: NUMERIC) -> List[int]:
        """Returns the indices of `data` equal to `compared`"""

        return np.where(data == compared)[0]

    def closest_relatives(strain: str, calls: pd.DataFrame,
                          distances: np.matrix) -> List[int]:
        """Return the row indices of the closest relatives of `strain`."""

        strain_index = tuple(calls.index).index(strain)
        strain_distances = np.array(distances[strain_index].flat)

        not_self = [strain_distances[:strain_index],
                    strain_distances[strain_index + 1:]]

        non_self_distances = np.concatenate(not_self)

        minimum_distance = min(non_self_distances)

        closest_indices = which(strain_distances, minimum_distance)

        # If the minimum dist is 0, it will still match to self here,
        # so ensure the query index is not in the return value
        closest = sorted(set(closest_indices) - {strain_index})

        return closest


    def closest_relative_allele(gene: str, closest_relative_indices,
                                calls: pd.DataFrame) -> List[int]:
        """For each of the genome indices specified by
        `closest_relative_indices`, return the allele found at `gene`.
        """
        closest_relative_alleles = calls[gene].iloc[closest_relative_indices]

        return closest_relative_alleles


    closest_indices = closest_relatives(strain, calls, distances)

    closest_alleles = closest_relative_allele(gene, closest_indices, calls)

    similarity = percent_shared(calls.loc[strain],
                                calls.iloc[closest_indices[0]])

    neighbour = Neighbour(closest_indices, closest_alleles, similarity)

    return neighbour


def find_locus_in_reference(gene: Path, reference: Path):

    def find_locus_location(query: str) -> int:
        """Executes a BLASTn search for a core gene against a reference genome.

        Returns the minimum of the start and stop locations of the gene
        so that gene calls can be ordered relative to the reference genome.
        """

        blast = ('blastn', '-task', 'megablast', '-subject', str(reference),
                 '-outfmt', '10')

        string_result = subprocess.run(blast, universal_newlines=True,
                                       check=True, input=query,
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

#@logtime('Reordering calls by reference')
def order_on_reference(reference: Path, genes: Path, calls: pd.DataFrame,
                       cores: int) -> pd.DataFrame:
    """Reorders `calls` columns to reflect the gene order found in `reference`.

    This is to facilitate analysis of gene linkage.
    """

    fasta_files = (gene for gene in genes.glob('*.fasta')
                   if gene.stem in calls.columns)

    with ProcessPoolExecutor(max_workers=cores) as ppe:

        find_locus = partial(find_locus_in_reference, reference=reference)

        gene_locations = list(ppe.map(find_locus, fasta_files))

    ordered_gene_locations = sorted(gene_locations, key=operator.itemgetter(1))

    ordered_gene_names, ordered_locs = zip(*ordered_gene_locations)

    return calls.reindex_axis(ordered_gene_names, axis=1)

#@logtime('Linkage disequilibrium')
def flank_linkage(strain: str, gene: str, hypothesis: int, gene_abundances,
                  calls: pd.DataFrame) -> Tuple[float, float]:
    """For three genes, (Left, Centre, Right), count how many observations of
    Centre were associated with the flanking genes Left and Right.
    """


    possible = np.array([a for a in gene_abundances[gene].keys()
                         if a != '?'])


    columns = tuple(calls.columns)
    gene_loc = columns.index(gene)

    left_col = columns[gene_loc - 1]

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

            flanks_given_h = (has_flank & is_h).sum() / has_flank.sum()

            if math.isnan(flanks_given_h):
                raise TypeError

        except TypeError:

            flanks_given_h = 0

    # If hypothesis is *never* observed with these flanks, it's not
    # impossible - merely unlikely. Use the probability of an unobserved allele
    # as the base probability
    return flanks_given_h or gene_abundances[gene]['?']


#@logtime('Partial sequence matching')
def partial_sequence_match(strain: str, gene: str, genes: Path,
                           jsondir: Path) -> Set[int]:
    """Attempts to use partial sequence data to exclude possible alleles."""

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


    def fragment_match(gene: str, fragment: str, genes: Path) -> Set[int]:
        """Attemps to match partial sequence data to a known allele from
        a multifasta file. Matches are only attemped at the beginning and end
        of the gene.
        """

        fragment_pattern = re.compile('(^{seq})|({seq}$)'.format(seq=fragment))

        glob_pattern = '*{}.f*'.format(gene)
        gene_file, *_ = genes.glob(glob_pattern)

        with gene_file.open('r') as fasta:

            result = set(int(rec.id) for rec in SeqIO.parse(fasta, 'fasta')
                         if re.search(fragment_pattern, str(rec.seq)))

        return result

    fragment = load_fragment(strain, gene, jsondir)

    matches = fragment_match(gene, fragment, genes)

    return matches

def allele_abundances(gene: str, calls: pd.DataFrame, replicates: int = 1000,
                      seed: int = 1) -> Dict[Union[str, int], float]:
    """Calculates the abundances of alleles for a given gene and attempts to
    determine the probability that the next observation will be a new allele.
    """

    observed_alleles = richness_estimate.Population(calls[gene])

    n_alleles = len(observed_alleles.abundance)

    discoveries = observed_alleles.monte_carlo(replicates, seed)

    last_percentile = max(1, int(0.01 * len(observed_alleles)))

    last_percentile_discovery_rate = mean(discoveries[-last_percentile:])

    to_subtract_from_known = last_percentile_discovery_rate / n_alleles

    abundance = {k: ((v / len(observed_alleles)) - to_subtract_from_known)
                 for k, v in observed_alleles.abundance.items()}

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

    adjusted_similarities = redistribute_allele_probability(combined_probs,
                                                            set(abundances))
    return adjusted_similarities

def bayes(adj_abundances, neighbour_probs, flanks) -> Dict[int, float]:


    def bayes_theorem(h: Union[str, int]) -> float:
        """Implementation of Bayes' Theorem in which the hypothesis being
        tested (h) is a given allele, and the lines of evidence are
        linkage disequilibrium between genes, and the similarity of the strain
        in question to its nearest neighbour.
        """

        neighbour_probability = neighbour_probs[h]

        flanks_given_h = flanks[h]

        p_h = adj_abundances[h]

        e_h = flanks_given_h * neighbour_probability

        return e_h * p_h

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

    neighbour = nearest_neighbour(gene, strain, fragment_matches,
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
                user_msg('Working on', gene, strain, '::', calls.loc[strain, gene])
                probs = recover_allele(strain, gene, calls, distances,
                                       genes, jsondir, gene_abundances)

                try:
                    results[strain][gene] = probs

                except KeyError:

                    results[strain] = {}
                    results[strain][gene] = probs

                counter += 1

                user_msg(counter, '/', total_bad)


def main():
    """Main function. Gathers arguments and passes them to recover()"""

    args = arguments()

    # diag
    recover(args.calls, args.reference, args.genes, args.jsons, args.distances,
            args.replicates, args.seed, args.cores)

if __name__ == '__main__':
    main()
