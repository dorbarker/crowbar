import argparse
import sys
from pathlib import Path

from .simulate_recovery import simulate_recovery
from .crossvalidate import crossvalidate


def arguments():

    parser = argparse.ArgumentParser()

    parser.set_defaults(func=None)
    subparsers = parser.add_subparsers(title='Commands')

    simulate_help = '''Simulate recovery by adding synthetic errors to a set of
                       known-good genomes.'''

    simulate = subparsers.add_parser('simulate',
                                     help=simulate_help)

    simulate.set_defaults(func=single_simulation)

    simulate.add_argument('--truncation-probability',
                        type=float,
                        default=0.0,
                        dest='trunc_prob',
                        help='Uniform probability that any given \
                              locus will be truncated [0.0]')

    simulate.add_argument('--missing-probability',
                        type=float,
                        default=0.0,
                        dest='miss_prob',
                        help='Uniform probability that any given \
                              locus will be rendered missing [0.0]')

    simulate.add_argument('--test-jsons',
                        type=Path,
                        required=True,
                        help='Directory containing FSAC-format JSONs')

    simulate.add_argument('--outdir',
                        type=Path,
                        required=True,
                        help='Output path')

    simulate.add_argument('--model',
                        type=Path,
                        required=True,
                        help='Path to pre-trained model')

    simulate.add_argument('-j', '--cores',
                        type=int,
                        default=1,
                        help='Number of CPU cores to use [1]')


    xval_help = '''Perform a cross-validation by simulating recovery in
                   replicate. A single pool of known-good genomes are used and
                   will be automatically divided into test and training sets.'''

    xval = subparsers.add_parser('crossvalidate', help=xval_help)

    xval.set_defaults(func=crossval)

    xval.add_argument('-n',
                      type=int,
                      required=True,
                      help='''The number of equally-sized subsamples into
                              which the population will be divided. For n
                              subsamples, n-1 will be used as training data.
                              ''')
    xval.add_argument('--seed',
                      type=int,
                      required=False,
                      default=None,
                      help='Seed to initialize the PRNG')

    xval.add_argument('--truncation-probability',
                      type=float,
                      default=0.0,
                      dest='trunc_prob',
                      help='Uniform probability that any given \
                            locus will be truncated [0.0]')

    xval.add_argument('--missing-probability',
                      type=float,
                      default=0.0,
                      dest='miss_prob',
                      help='Uniform probability that any given \
                            locus will be rendered missing [0.0]')

    xval.add_argument('--jsons',
                      type=Path,
                      required=True,
                      help='Directory containing FSAC-formatted JSONs')

    xval.add_argument('--outdir',
                      type=Path,
                      required=True,
                      help='Output directory')

    xval.add_argument('-a', '--alleles',
                      type=PAth,
                      required=True,
                      help='Path to directory containing FASTA-format alleles')

    xval.add_argument('-j', '--cores',
                      type=int,
                      required=False,
                      help='Number of CPU cores to use [1]')

    args = parser.parse_args()
    if args.func is None:
        parser.print_help()
        sys.exit(0)

    return args


def single_simulation(args):

    simulate_recovery(args.outdir, args.test_jsons, args.model,
                      args.trunc_prob, args.miss_prob, args.cores)


def crossval(args):

    crossvalidate.crossvalidate(args.n, args.trunc_prob, args.miss_prob,
                                args.outdir, args.jsons, args.alleles,
                                args.cores, args.seed)

def main():

    args = arguments()

    args.func(args)

if __name__ == '__main__':
    main()

