import argparse
import shutil

from pathlib import Path

import pandas as pd


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
    ...


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


def calculate_abundances():
    ...


def save_abundances():
    ...


def calculate_distance_matrix():
    ...


def save_distance_matrix():
    ...


if __name__ == '__main__':
    main()
