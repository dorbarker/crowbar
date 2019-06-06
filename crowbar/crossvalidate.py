import math
import itertools
import random
from pathlib import Path
from typing import Optional

from simulate_recovery import PathTable

import model

print(model.build_model)

def crossvalidate():
    pass


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


def generate_models(paths: PathTable):

    for experiment in paths['experiments'].glob('*/'):
        pass


def run_experiments():
    pass

