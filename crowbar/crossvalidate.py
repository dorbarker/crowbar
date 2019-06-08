import math
import itertools
import random
from pathlib import Path
from typing import Optional
import subprocess

from .simulate_recovery import PathTable

from . import model
from . import simulate_recovery


def crossvalidate(n: int, trunc_prob: float, miss_prob: float,
                  experiments: Path, json_dir: Path, alleles: Path, cores: int,
                  seed: Optional[int] = None):

    paths = generate_paths_table(experiments, json_dir, alleles)

    divide_jsons(n, paths, seed)

    generate_models(paths, cores)

    run_experiments(paths, trunc_prob, miss_prob, cores)


def generate_paths_table(experiments: Path, jsons: Path,
                         alleles: Path) -> PathTable:

    paths = {
        'experiments': experiments,
        'alleles':     alleles,
        'jsons':       jsons
    }

    paths['experiment'].mkdir(parents=True)

    return paths


def divide_jsons(n: int, paths: PathTable, seed: Optional[int] = None) -> None:

    random.seed(seed)

    jsons = list(paths['jsons'].glob('*.json'))

    random.shuffle(jsons)

    chunk_length = math.ceil(len(jsons) / n)

    json_range = range(0, len(jsons), chunk_length)

    experiment_groups = [jsons[i:i+chunk_length] for i in json_range]

    experiments = itertools.product(range(n),  enumerate(experiment_groups))

    for experiment, (group_index, json_group) in experiments:
        # For each json group, create a directory where it is the test set, and
        # add them to the training set of all other groups

        d = paths['experiments'].joinpath(f'experiment_{experiment}')

        test = d.joinpath('test')
        train = d.joinpath('training')

        if experiment == group_index:

            dst = test
        else:

            dst = train

        dst.mkdir(parents=True, exist_ok=True)

        for j in json_group:

            out_json = dst.joinpath(j.name)

            out_json.symlink_to(j.resolve())


def generate_models(paths: PathTable, cores: int):

    for experiment in paths['experiments'].glob('*/'):

        training_calls = experiment / 'training_calls'

        # TODO consider changing subprocess.run to a proper fsac module import
        tabulate_training_calls = ('fsac', 'tabulate',
                                   '--json-dir',  experiment / 'training',
                                   '--output',    training_calls,
                                   '--delimiter', ',')

        subprocess.run(tabulate_training_calls, check=True)

        calls = model.load_calls(training_calls)

        model.build_model(calls, paths['alleles'], experiment / 'model', cores)


def run_experiments(paths: PathTable, trunc_prob, miss_prob, cores):

    for experiment in paths['experiments'].glob('*/'):

        simulate_recovery.simulate_recovery(experiment / 'results',
                                            experiment / 'test',
                                            experiment / 'model',
                                            trunc_prob,
                                            miss_prob,
                                            cores)
