#!/usr/bin/env python3

import argparse
import json
import sys

from collections import namedtuple
from operator import itemgetter
from pathlib import Path

import pandas as pd


AlleleProb = namedtuple('AlleleProb', ('allele', 'probability'))


def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('--recovered',
                        help='crowBAR output')

    parser.add_argument('--calls',
                        type=Path,
                        required=True,
                        help='Table of allele calls')

    parser.add_argument('-o', '--output',
                        type=Path,
                        help='Updated calls table')

    return parser.parse_args()


def main():

    args = arguments()

    data = load_data(args.recovered)

    results_table = tabulate_results(data)

    updated_calls = update_calls(data, args.calls)

    results_table.to_csv(args.output.parent / 'probability_ratios.csv')

    updated_calls.to_csv(args.output)


def load_data(inpath):

    file_obj = inpath or sys.stdin

    try:
        with file_obj.open('r') as f:
            data = json.load(f)

    except AttributeError:
        data = json.load(file_obj)

    finally:
        return data


def iter_data(data):

    for strain in data:
        for gene in data[strain]:

            yield strain, gene, data[strain][gene]


def tabulate_results(data) -> pd.DataFrame:

    def best_two(recovery):

        ordered = sorted(recovery.items(), key=itemgetter(1), reverse=True)

        try:

            first, second, *_ = ordered

        except ValueError:

            first, second = [*ordered, (None, 0)]

        return AlleleProb(*first), AlleleProb(*second)

    rows = []

    for strain, gene, value in iter_data(data):

        fst, snd = best_two(value)

        try:
            ratio = fst.probability / snd.probability

        except ZeroDivisionError:

            ratio = None

        row = {
            'strain': strain,
            'gene': gene,
            'best': fst.allele,
            'best_probability': fst.probability,
            'second': snd.allele,
            'second_probability': snd.probability,
            'ratio': ratio
        }

        rows.append(row)

    table = pd.DataFrame(rows)

    headers = ('strain', 'gene', 'best', 'best_probability',
               'second', 'second_probability', 'ratio')

    return table.reindex_axis(headers, axis=1)


def update_calls(data, calls_path: Path) -> pd.DataFrame:

    calls = pd.read_csv(calls_path)

    for strain, gene, probabilities in iter_data(data):

        most_probable = max(probabilities, key=lambda x: probabilities[x])

        calls.loc[strain, gene] = most_probable

    return calls


if __name__ == '__main__':
    main()
