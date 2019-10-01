import argparse
import logging
import sys
from pathlib import Path

from . import __version__
from . import verbose_help

from .recover import recover_genome
from .model import load_calls, build_model

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('-v', '--version',
                        action='version',
                        version=f'{parser.prog} {__version__}')

    parser.set_defaults(func=None)
    subparsers = parser.add_subparsers(title='Commands')

    ### Recovery
    recovery_help = 'Recover untypable alleles'

    recover = subparsers.add_parser('recover', help=recovery_help)
    recover.set_defaults(func=recover_alleles)

    recover.add_argument('-i', '--input',
                         type=Path,
                         required=True,
                         help='Path to input fsac-formatted JSON')

    recover.add_argument('-o', '--output',
                         type=Path,
                         required=True,
                         help='Path to output directory')

    recover.add_argument('--model',
                         type=Path,
                         required=True,
                         help='Path to recovery model')

    recover.add_argument('-V', '--verbose',
                         action='store_true',
                         help=verbose_help)

    ### Model Building
    model_help = 'Build a recovery model from allele call data'

    build_model = subparsers.add_parser('model', help=model_help)
    build_model.set_defaults(func=train_model)

    build_model.add_argument('--calls',
                             type=Path,
                             required=True,
                             help='Path to allele calls table')

    build_model.add_argument('--alleles-dir',
                             type=Path,
                             required=True,
                             help='Path to directory containing gene FASTAs')

    build_model.add_argument('--output',
                             type=Path,
                             required=True,
                             help='Directory containing model')

    build_model.add_argument('--cores',
                             type=int,
                             default=1,
                             help='Number of CPU cores to use [1]')

    build_model.add_argument('-V', '--verbose',
                             action='store_true',
                             help=verbose_help)


    args = parser.parse_args()

    if args.func is None:
        parser.print_help()
        sys.exit(0)

    return args


def recover_alleles(args):

    recover_genome(args.input, args.model, args.output)


def train_model(args):

    calls = load_calls(args.calls)

    build_model(calls, args.alleles_dir, args.output, args.cores)


def main():

    args = arguments()

    logging.basicConfig(datefmt='%Y-%m-%d %H:%M',
                        format='%(asctime)s - %(levelname)s: %(message)s',
                        stream=sys.stderr,
                        level=logging.INFO if args.verbose else logging.ERROR)

    logging.info("Running %s", args.func.__name__)

    args.func(args)

    logging.info("%s complete", args.func.__name__)

if __name__ == '__main__':
    main()

