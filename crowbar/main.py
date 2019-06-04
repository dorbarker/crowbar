import argparse
import sys

from . import __version__

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('-v', '--version',
                        action='store_true',
                        help='Print version and exit')

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


    args = parser.parse_args()

    if args.version:
        print('crowbar', __version__)
        sys.exit(0)

    if args.func is None:
        parser.print_help()
        sys.exit(0)

    return args


def main():

    args = arguments()

    args.func(args)


def recover_alleles():
    ...


def train_model():
    ...

if __name__ == '__main__':
    main()

