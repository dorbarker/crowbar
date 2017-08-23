import argparse
import json
import re
import operator
import subprocess
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from functools import reduce
from io import StringIO
from pathlib import Path
from statistics import mean
from typing import Callable, Dict, List, Sequence, Set, Tuple, Union

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

    parser = argparse.ArgumentParser()

    parser.add_argument('--cores',
                        type=int,
                        default=1,
                        help='Number of CPU cores to use [1]')

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

    def non_missing_hamming(j):

        strain1 = row[1]
        strain2 = calls.iloc[j]

        shared = [a and b for a, b in zip(strain1 > 0, strain2 > 0)]

        return sum(strain1[shared] != strain2[shared])

    return {j: non_missing_hamming(j) for j in range(idx + 1, len(calls))}


def dist_gene(calls: pd.DataFrame, cores: int) -> np.matrix:


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


def closest_relatives(strain: str, calls: pd.DataFrame,
                      distances: np.matrix) -> List[int]:

    def which(data: Sequence, operator_: Callable,
              compared: NUMERIC) -> List[int]:

        return [i for i, v in enumerate(data) if operator_(v, compared)]

    query_index = calls.columns.index(strain)
    query_distances = distances[query_index]

    non_self_distances = query_distances[:query_index] + \
                         query_distances[query_index + 1:]

    minimum_distance = min(non_self_distances)

    closest_indices = which(query_distances, operator.eq, minimum_distance)

    # If the minimum dist is 0, it will still match to self here,
    # so ensure the query index is not in the return value
    closest = sorted(set(closest_indices) - {query_index})

    return closest


def closest_relative_allele(gene: str, closest_relative_indices,
                            calls: pd.DataFrame) -> List[int]:

    closest_relative_alleles = calls[gene].iloc[closest_relative_indices]

    return closest_relative_alleles


def percent_shared(strain1, strain2):

    shared = [a and b for a, b in zip(strain1 > 0, strain2 > 0)]

    return sum(strain1[shared] == strain2[shared]) / len(shared)


def exclude_matches(include: Set[int], matches: Sequence[int]) -> List[int]:

    return [match for match in matches if match in include]


def order_on_reference(reference: Path, genes: Path,
                       calls: pd.DataFrame) -> pd.DataFrame:

    gene_locations = []

    for gene in genes.glob('*.fasta'):

        gene_name = gene.stem

        with gene.open('r') as fasta:

            # Use just the first record
            rec = next(SeqIO.parse(fasta, 'fasta'))

            query = str(rec.seq)

            loc = find_locus_location(query, reference)

            gene_locations.append((gene_name, loc))

    ordered_gene_locations = sorted(gene_locations, key=operator.itemgetter(1))

    ordered_gene_names, ordered_locs = zip(*ordered_gene_locations)

    return calls.reindex_axis(ordered_gene_names, axis=1)


def find_locus_location(query: str, subject: Path) -> int:

    blast = ('blastn', '-task', 'megablast', '-subject', str(subject),
             '-outfmt', '10')

    string_result = subprocess.run(blast, universal_newlines=True, check=True,
                                   input=query,
                                   stdout=subprocess.PIPE)

    table_result = pd.DataFrame(StringIO(string_result.stdout))

    # best row should be first
    # sstart, ssend = 8, 9
    start, stop = table_result.iloc[0, [8, 9]]

    return min(start, stop)


def count_triplets(gene: str, calls: pd.DataFrame) -> TRIPLET_COUNTS:

    def tree():
        return defaultdict(tree)

    left_col = calls.columns[calls.columns.index(gene) - 1]
    right_col = calls.columns[calls.columns.index(gene) + 1]

    triplet_calls = calls[[left_col, gene, right_col]]

    triplets = tree()

    for _, row in triplet_calls.iterrows():

        if -1 in row or 0 in row:
            continue


        left, centre, right = row

        try:
            triplets[left][right][centre] += 1

        except TypeError:  # if not yet initialized
            triplets[left][right][centre] = 1

    # Convert back to standard dict to avoid any weirdness later
    return json.loads(json.dumps(triplets))


def fragment_match(gene: str, fragment: str, genes: Path) -> Set[int]:

    fragment_pattern = re.compile('^{seq}|{seq}$'.format(seq=fragment))

    glob_pattern = '*{}.f*'.format(gene)
    gene_file, *_ = genes.glob(glob_pattern)

    with gene_file.open('r') as fasta:

        result = set(rec.id for rec in SeqIO.parse(fasta, 'fasta')
                     if re.search(fragment_pattern, rec.seq))

    return result

def load_fragment(strain: str, gene: str, jsondir: Path) -> str:
    # Currently, totally MIST-based

    jsonpath = (jsondir / strain).with_suffix('.json')

    with jsonpath.open('r') as json_obj:
        data = json.load(json_obj)

    test_results = data['Results'][0]['TestResults']

    # presumed to be a single test name per file
    test_name, *_ = test_results.keys()

    test_data = test_results[test_name]

    fragment = test_data[gene]['Amplicon']

    return fragment


def linkage_disequilibrium(gene: str, strain: str,
                           calls: pd.DataFrame) -> TRIPLET_COUNT_CURRENT:

    gene_names = calls.columns

    triplets = count_triplets(gene, calls)

    left_col = gene_names[gene_names.index(gene) - 1]
    right_col = gene_names[gene_names.index(gene) + 1]

    left, right = calls[[left_col, right_col]].loc[strain]

    strain_flanking = triplets[left][right]

    return triplets, strain_flanking


def linkage_probability(triplets, fragment_matches) -> TRIPLET_PROBS:

    filtered_linkage = {k: v for k, v in triplets.items()}

    linkage_probabilities = {}  # type: Dict[int, Dict]

    total = 0

    for left in filtered_linkage:
        for right in filtered_linkage[left]:
            for centre in filtered_linkage[left][right]:

                if centre not in fragment_matches:
                    del filtered_linkage[left][right][centre]
                else:
                    total += triplets[left][right][centre]

    for left in filtered_linkage:
        linkage_probabilities[left] = {}
        for right in filtered_linkage[left]:
            linkage_probabilities[left][right] = {}

            for centre in filtered_linkage[left][right]:

                count = filtered_linkage[left][right][centre]

                linkage_probabilities[left][right][centre] = count / total

    return linkage_probabilities


def nearest_neighbour(gene: str, strain: str, included_fragments: Set[int],
                      distances: np.matrix, calls: pd.DataFrame) -> List[int]:

    closest_alleles = closest_relative_allele(strain, gene, calls, distances)

    filtered_matches = exclude_matches(included_fragments, closest_alleles)

    return filtered_matches


def partial_sequence_match(gene: str, strain: str, genes: Path,
                           jsondir: Path) -> Set[int]:

    fragment = load_fragment(gene, strain, jsondir)

    matches = fragment_match(gene, fragment, genes)

    return matches


def allele_abundances(gene: str, calls: pd.DataFrame, replicates: int = 1000,
                      seed: int = 1) -> Dict[Union[str, int], float]:

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
        return previous + (current * (1 - previous))

    allele_proportions = Counter(neighbour_alleles)

    combined_probs = {k: reduce(combine, (abundances[k] for _ in range(count)))
                      for allele, count in allele_proportions.items()}

    inverses = {k: (1 - neighbour_similarity) * (1 - abundances[k])
                for k in combined_probs}

    return combined_probs, inverses

def bayes(abundances: Dict[int, float], fragment_matches: Set[int],
          triplets: TRIPLET_COUNTS, linkage: Dict[int, int],
          neighbour_similarity: float,
          neighbour_alleles: List[int]) -> Dict[int, float]:


    # the priors
    adj_abundances = redistribute_next_allele_probability(abundances,
                                                          fragment_matches)

    # linkage counts coverted to probabilities
    linkage_probs = linkage_probability(triplets, fragment_matches)

    probs_inverse = combine_neighbour_similarities(neighbour_similarity,
                                                   neighbour_alleles,
                                                   adj_abundances)

    neighbour_probs, inverse_neighbour = probs_inverse

    def marginal_likelihood(h):

        def neighbour_terms():

            terms = [neighbour_probs[k] * inverse_neighbour[k]
                     for k in neighbour_probs]

            return terms

        linkage_total = sum(linkage.values())

        not_h_abundance = 1 - adj_abundances[h]

        h_with_flanks = linkage_probs[h]
        not_h_with_flanks = 1 - h_with_flanks

        linkage_likelihood = \
                ((h_with_flanks / linkage_total) * adj_abundances[h]) + \
                ((not_h_with_flanks / linkage_total) * not_h_abundance)


        return reduce(operator.mul, [linkage_likelihood] + neighbour_terms)


    def bayes_theorem(h):

        if h in neighbour_alleles:

            neighbour = neighbour_probs[h]

        else:
            # neighbour = 1 - (neighbour_similarity * adj_abundances[h])
            neighbour = inverse_neighbour[h]

        p_h = adj_abundances[h]

        e = marginal_likelihood(h)

        e_h = (triplets[h] / sum(triplets.values())) * neighbour

        h_e = (e_h * adj_abundances[h]) / e

    return {h: bayes_theorem(h) for h in adj_abundances}

def recover_allele(strain: str, gene: str, calls: pd.DataFrame,
                   distances: np.matrix, genes: Path, jsondir: Path,
                   replicates: int, seed: int):

    if calls.loc[strain, gene] == -1:  # is truncation

        fragment = load_fragment(strain, gene, jsondir)

        fragment_matches = fragment_match(gene, fragment, genes)

    else:  # is wholly missing

        fragment_matches = set(calls[gene])

    nearest_neighbours = nearest_neighbour(gene, strain, fragment_matches,
                                           distances, calls)

    neighbour_alleles = closest_relative_allele(gene, nearest_neighbours,
                                                calls)

    neighbour_similarity = percent_shared(calls.loc[strain],
                                          calls.iloc[nearest_neighbours[0]])

    triplets, strain_flanking = linkage_disequilibrium(gene, strain, calls)

    abundance = allele_abundances(gene, calls, replicates, seed)

    probs = bayes(abundance, fragment_matches, triplets, strain_flanking,
                  neighbour_similarity, neighbour_alleles)

    return probs

def recover(callspath: Path, reference: Path, genes: Path, jsondir: Path,
            replicates: int, seed: int, cores: int):

    calls = order_on_reference(reference, genes,
                               calls=pd.read_csv(callspath, index_col=0))

    distances = dist_gene(calls, cores)

    for strain in calls.index:
        for gene in calls.columns:

            if not calls.loc[strain, gene] > 0:  # truncated or missing
                recover_allele(strain, gene, calls, distances, genes, jsondir,
                               replicates, seed)


def main():

    args = arguments()

if __name__ == '__main__':
    main()
