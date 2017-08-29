import argparse
import json
import re
import operator
import subprocess
from collections import Counter, defaultdict, namedtuple
from concurrent.futures import ProcessPoolExecutor
from functools import reduce
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


def row_distance(idx: int, row: Tuple[str, pd.Series],
                 calls: pd.DataFrame) -> Dict[int, int]:
    """Returns the Hamming distance of non-missing alleles between two strains.

    Results are returned as a dictionary for distances between the query strain
    (strain1) and all the subject strain (each strain2).

    Called by dist_gene()
    """

    def non_missing_hamming(j):
        """Returns the distance between two strains, considering only loci
        which are not missing in either individual.
        """

        strain1 = row[1]
        strain2 = calls.iloc[j]

        shared = [a and b for a, b in zip(strain1 > 0, strain2 > 0)]

        return sum(strain1[shared] != strain2[shared])

    return {j: non_missing_hamming(j) for j in range(idx + 1, len(calls))}


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

        dist_mat = np.matrix([[0 for _ in range(n_row)] for _ in range(n_row)])

        with ProcessPoolExecutor(max_workers=cores) as ppe:
            futures = {i: ppe.submit(row_distance, i, row, calls)
                       for i, row in enumerate(calls.iterrows())}

        results = {i: j.result() for i, j in futures.items()}

        for i, js in results.items():

            for j in js:
                dist_mat[i, j] = dist_mat[j, i] = results[i][j]

        return dist_mat


    try:

        distances = np.matrix(pd.read_csv(distance_path))

    except FileNotFoundError:

        distances = dist_gene(calls, cores)
        pd.DataFrame(distances).to_csv(distance_path)

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


    def which(data: Sequence, operator_: Callable,
              compared: NUMERIC) -> List[int]:
        """Returns the indices of `data` for which `operator_` is True"""

        return [i for i, v in enumerate(data) if operator_(v, compared)]

    def closest_relatives(strain: str, calls: pd.DataFrame,
                          distances: np.matrix) -> List[int]:
        """Return the row indices of the closest relatives of `strain`."""

        strain_index = tuple(calls.index).index(strain)
        strain_distances = tuple(distances[strain_index].flat)

        non_self_distances = strain_distances[:strain_index] + \
                             strain_distances[strain_index + 1:]

        minimum_distance = min(non_self_distances)

        closest_indices = which(strain_distances,
                                operator.eq, minimum_distance)

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

    filtered_matches = exclude_matches(included_fragments, closest_alleles)

    similarity = percent_shared(calls.loc[strain],
                                calls.iloc[closest_indices[0]])

    neighbour = Neighbour(closest_indices, closest_alleles, similarity)

    return neighbour


def exclude_matches(include: Set[int], matches: Sequence[int]) -> List[int]:
    """Filters `matches` based on their presence in `include`"""

    return [match for match in matches if match in include]


def order_on_reference(reference: Path, genes: Path,
                       calls: pd.DataFrame) -> pd.DataFrame:
    """Reorders `calls` columns to reflect the gene order found in `reference`.

    This is to facilitate analysis of gene linkage.
    """

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


    gene_locations = []

    for gene in genes.glob('*.fasta'):

        gene_name = gene.stem

        with gene.open('r') as fasta:

            # Use just the first record
            rec = next(SeqIO.parse(fasta, 'fasta'))

            query = str(rec.seq)

            loc = find_locus_location(query)

            gene_locations.append((gene_name, loc))

    ordered_gene_locations = sorted(gene_locations, key=operator.itemgetter(1))

    ordered_gene_names, ordered_locs = zip(*ordered_gene_locations)

    return calls.reindex_axis(ordered_gene_names, axis=1)


def count_triplets(gene: str, calls: pd.DataFrame) -> TRIPLET_COUNTS:
    """For three genes, (Left, Centre, Right), count how many observations of
    Centre were associated with the flanking genes Left and Right.
    """

    def tree():
        """Factory function for generating tree-like dict structures"""

        return defaultdict(tree)

    incomplete = pd.Series([0, -1])

    columns = tuple(calls.columns)
    gene_loc = columns.index(gene)

    left_col = columns[gene_loc - 1]
    right_col = columns[gene_loc + 1]

    triplet_calls = calls[[left_col, gene, right_col]]

    triplets = tree()

    for _, row in triplet_calls.iterrows():

        if any(incomplete.isin(row)):
            continue

        left, centre, right = row

        try:
            triplets[left][right][centre] += 1

        except TypeError:  # if not yet initialized
            triplets[left][right][centre] = 1

    return triplets

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

        fragment_pattern = re.compile('^{seq}|{seq}$'.format(seq=fragment))

        glob_pattern = '*{}.f*'.format(gene)
        gene_file, *_ = genes.glob(glob_pattern)

        with gene_file.open('r') as fasta:

            result = set(int(rec.id) for rec in SeqIO.parse(fasta, 'fasta')
                         if re.search(fragment_pattern, str(rec.seq)))

        return result

    fragment = load_fragment(strain, gene, jsondir)

    matches = fragment_match(gene, fragment, genes)

    return matches


def linkage_disequilibrium(gene: str, strain: str,
                           calls: pd.DataFrame) -> TRIPLET_COUNT_CURRENT:

    gene_names = tuple(calls.columns)

    triplets = count_triplets(gene, calls)

    left_col = gene_names[gene_names.index(gene) - 1]
    right_col = gene_names[gene_names.index(gene) + 1]

    left, right = calls[[left_col, right_col]].loc[strain]

    strain_flanking = triplets[left][right]

    return triplets, strain_flanking


#def linkage_probability(triplets, fragment_matches) -> TRIPLET_PROBS:
#
#    filtered_linkage = {k: v for k, v in triplets.items()}
#
#    linkage_probabilities = {}  # type: Dict[int, Dict]
#
#    total = 0
#
#    for left in filtered_linkage:
#        for right in filtered_linkage[left]:
#
#            for centre in filtered_linkage[left][right]:
#
#                to_del = set()
#
#                if centre not in fragment_matches:
#                    to_del.add(centre)
#                else:
#                    total += triplets[left][right][centre]
#
#            l_r_items = filtered_linkage[left][right].items()
#
#            new_centre = {k: v for k, v in l_r_items if k not in to_del}
#
#            filtered_linkage[left][right] = new_centre
#
#    for left in filtered_linkage:
#
#        for right in filtered_linkage[left]:
#
#            right_candidates = {}
#
#            for centre in filtered_linkage[left][right]:
#
#                count = filtered_linkage[left][right][centre]
#
#                try:
#                    right_candidates[centre] = count / total
#
#                except ZeroDivisionError:
#                    right_candidates[centre] = 0.0
#
#            if right_candidates:
#
#                try:
#
#                    linkage_probabilities[left][right] = right_candidates
#
#                except KeyError:
#
#                    linkage_probabilities[left] = {}
#
#                    linkage_probabilities[left][right] = right_candidates
#
#    return linkage_probabilities
#
def linkage_probability(flanks: Dict[int, int],
                        fragment_matches: Set[int]) -> Dict[int, float]:

    filtered_flanks = {k: v for k, v in flanks.items()
                       if k in fragment_matches}

    sum_remaining = sum(filtered_flanks.values())

    linkage_probabilities = {k: v / sum_remaining
                             for k, v in filtered_flanks.items()}

    return linkage_probabilities


def allele_abundances(gene: str, calls: pd.DataFrame, replicates: int = 1000,
                      seed: int = 1) -> Dict[Union[str, int], float]:
    """Calculates the abundances of alleles for a given gene and attempts to
    determine the probability that the next observation will be a new allele.
    """

    observed_alleles = richness_estimate.Population(calls[gene])

    n_alleles = len(observed_alleles.abundance)

    discoveries = observed_alleles.monte_carlo(replicates, seed)

    last_percentile = int(0.01 * len(observed_alleles))

    last_percentile_discovery_rate = mean(discoveries[-last_percentile:])

    to_subtract_from_known = last_percentile_discovery_rate / n_alleles

    abundance = {k: ((v / len(observed_alleles)) - to_subtract_from_known)
                 for k, v in observed_alleles.abundance.items()}

    abundance['?'] = last_percentile_discovery_rate

    return abundance


def redistribute_next_allele_probability(abundances, fragment_matches):
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


def combine_neighbour_similarities(neighbour_similarity: float,
                                   neighbour_alleles: List[int],
                                   abundances: Dict[Union[int, str], float]):

    def combine(previous, current):
        """Multiplies the inverse probability of the previous iteration by
        the probability of the current iteration.
        """

        return previous + (current * (1 - previous))

    allele_proportions = Counter(neighbour_alleles)
    print('Proportions:', allele_proportions)
    combined_probs, inverses = {}, {}

    for k, count in allele_proportions.items():

        try:
            obs = (abundances[k] for _ in range(count))
            combined_probs[k] = reduce(combine, obs)

            inverses[k] = (1 - neighbour_similarity) * (1 - abundances[k])

        except KeyError:
            combined_probs[k] = neighbour_similarity
            inverses[k] = 1.0 - neighbour_similarity

    return combined_probs, inverses


def bayes(abundances: Dict[Union[str, int], float], fragment_matches: Set[int],
          triplets: TRIPLET_COUNTS, linkage: Dict[int, int],
          neighbour: Neighbour) -> Dict[int, float]:


    # the priors
    adj_abundances = redistribute_next_allele_probability(abundances,
                                                          fragment_matches)

    # linkage counts coverted to probabilities
    linkage_probs = linkage_probability(linkage, fragment_matches)

    probs_inverse = combine_neighbour_similarities(neighbour.similarity,
                                                   neighbour.alleles,
                                                   adj_abundances)

    neighbour_probs, inverse_neighbour = probs_inverse

    def marginal_likelihood(h: Union[str, int]) -> float:

        def prepare_neighbour_terms() -> List[float]:
            """Returns a list of terms for the denominator of Bayes Theorem."""

            terms = [neighbour_probs[k] * inverse_neighbour[k]
                     for k in neighbour_probs]

            return terms

        not_h_abundance = 1 - adj_abundances[h]

        try:

            h_with_flanks = linkage_probs[h]

        except KeyError:

            h_with_flanks = 0.0

        not_h_with_flanks = 1 - h_with_flanks

        linkage_likelihood = \
                (h_with_flanks * adj_abundances[h]) + \
                (not_h_with_flanks * not_h_abundance)

        neighbour_terms = prepare_neighbour_terms()

        # ... I think?
        terms = [term for term in (linkage_likelihood, *neighbour_terms)
                 if term > 0]

        print('terms:', terms)
        return reduce(operator.mul, terms)


    def bayes_theorem(h: Union[str, int]) -> float:
        """Implementation of Bayes' Theorem in which the hypothesis being
        tested (h) is a given allele, and the lines of evidence are
        linkage disequilibrium between genes, and the similarity of the strain
        in question to its nearest neighbour.
        """
        print('Hypothesis:', h)
        print('Inverse:', inverse_neighbour)
        print('Neighbour:', neighbour_probs)
        if h in neighbour.alleles:
            neighbour_prob = neighbour_probs[h]

        # TODO: double check this
        else:

            neighbour_prob = 1.0 - neighbour.similarity

        p_h = adj_abundances[h]

        e = marginal_likelihood(h)

        try:

            e_h = (linkage[h] / \
                   sum(int(x) for x in linkage.values() if isinstance(x, int))) \
                   * neighbour_prob

        except TypeError:

            e_h = neighbour_prob

        h_e = (e_h * p_h) / e

        print('e_h:', e_h)
        print('p_h:', p_h)
        print('e:', e)

        print('Probability:', h_e, '\n\n')
        return h_e

    print('Adjusted abundances:', adj_abundances)
    return {h: bayes_theorem(h) for h in adj_abundances}


def recover_allele(strain: str, gene: str, calls: pd.DataFrame,
                   distances: np.matrix, genes: Path, jsondir: Path,
                   replicates: int, seed: int):
    """Attempts to recover the allele of a missing locus based on
    lines of evidence from allele abundance, linkage disequilibrium, the
    close relatives.
    """
    if calls.loc[strain, gene] == -1:  # is truncation

        fragment_matches = partial_sequence_match(strain, gene, genes, jsondir)

    else:  # is wholly missing

        fragment_matches = set(calls[gene])

    neighbour = nearest_neighbour(gene, strain, fragment_matches,
                                  distances, calls)

    abundance = allele_abundances(gene, calls, replicates, seed)

    triplets, strain_flanking = linkage_disequilibrium(gene, strain, calls)

    probs = bayes(abundance, fragment_matches, triplets, strain_flanking,
                  neighbour)

    return probs


def recover(callspath: Path, reference: Path, genes: Path, jsondir: Path,
            distance_path: Optional[Path], replicates: int, seed: int,
            cores: int):
    """Master function for crowBAR.

    Walks through a pandas DataFrame of allele calls and attempts to recover
    any truncated (-1) or absent (0) loci.
    """


    results = {}  # Dict[str, Dict[str, Dict[int, float]

    calls = order_on_reference(reference, genes,
                               calls=pd.read_csv(callspath, index_col=0))

    distances = hamming_distance_matrix(distance_path, calls, cores)
    print("Done getting order and distances")

    for strain in calls.index:
        results[strain] = {}
        for gene in calls.columns:

            if not calls.loc[strain, gene] > 0:  # truncated or missing
                print('Working on', gene, strain, '::', calls.loc[strain, gene])
                probs = recover_allele(strain, gene, calls, distances,
                                       genes, jsondir, replicates, seed)

                results[strain][gene] = probs


    print(results)


def main():

    args = arguments()

    # diag
    recover(args.calls, args.reference, args.genes, args.jsons, args.distances,
            1000, 1, args.cores)

if __name__ == '__main__':
    main()
