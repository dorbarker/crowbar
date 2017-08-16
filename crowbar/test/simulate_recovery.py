# import from one level above
import sys
import os.path

up = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
sys.path.append(up)

import crowbar  # main script

# regular imports
import argparse
import random
from pathlib import Path
from typing import Dict, Tuple
# 3rd Party imports
import pandas as pd
def arguments():

    parser = argparse.ArgumentParser()

    parser = parser.add_argument('--seed',
                                 type=int,
                                 default=1,
                                 help='Seed to initialize the PRNG [1]')

    parser = parser.add_argument('--truncation-probability',
                                 type=float,
                                 default=0.0,
                                 store='trunc_prob',
                                 help='Uniform probability that any given \
                                       locus will be truncated [0.0]')

    parser = parser.add_argument('--missing-probability',
                                 type=float,
                                 default=0.0,
                                 store='miss_prob',
                                 help='Uniform probability that any given \
                                       locus will be rendered missing [0.0]')

    return parser.parse_args()

def sequential_test():

    # for each calls[strain, gene], truncate or vanish
        # get results
        # return most probable match
    pass

def introduce_random_errors(trunc_prob: float, miss_prob: float,
                            calls: pd.DataFrame, jsondir: Path,
                            seed: int) -> Tuple[pd.DataFrame, Dict[str, str]]:

    random.seed(seed)

    for strain, row in calls.iterrows():

        to_truncate = [random.random() < trunc_prob for _ in row]
        to_vanish = [random.random() < miss_prob for _ in row]

        row[to_truncate] = [-1 for _ in row]
        row[to_vanish] = [0 for _ in row]

        truncations = {gene: truncate(strain, gene, jsondir)
                       for gene in row[to_truncate].index}

    return calls, truncations
def truncate(strain: str, gene: str, jsondir: Path):

    sequence = crowbar.load_fragment(strain, gene, jsondir)

    # Leave at least a 50-mer on each end
    pivot = random.randint(50, len(sequence) - 50)
    side = random.choice((0, 1))

    halves = sequence[:pivot], sequence[:pivot]

    return halves[side]

def vanish(strain: str, gene: str):
    pass

def main():

    args = arguments()

if __name__ == '__main__':
    main()
