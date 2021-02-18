# Copyright 2019-2020 John Lees

'''CLI for sparrowhawk'''

import os, sys

import sparrowhawk

from .__init__ import __version__

def get_options():
    import argparse

    description = 'Assembly microbial genomes using GPUs'
    parser = argparse.ArgumentParser(description=description,
                                     prog='sparrohawk')

    modeGroup = parser.add_argument_group('Mode of operation')
    mode = modeGroup.add_mutually_exclusive_group(required=True)
    mode.add_argument('--assemble',
                        action='store_true',
                        default=False,
                        help='Assemble a genome')
    mode.add_argument('--dbg',
                        action='store_true',
                        default=False,
                        help='Build a DBG at a single k')
    mode.add_argument('--info',
                        action='store_true',
                        default=False,
                        help='Return information about the system and the software')

    io = parser.add_argument_group('Input/output')
    io.add_argument('--fwd',
                    required=True,
                    help='Forward reads')
    io.add_argument('--rev',
                    required=True,
                    help='Reverse reads')
    io.add_argument('--output',
                    required=True,
                    default='sparrowhawk',
                    help="Output prefix [default = 'sparrowhawk']")

    kmerGroup = parser.add_argument_group('Kmer options')
    kmerGroup.add_argument('--k', required=True,
                            help='Comma separated k-mer lengths')

    optimisation = parser.add_argument_group('Optimisation options')
    optimisation.add_argument('--cpus',
                              type=int,
                              default=1,
                              help='Number of CPUs to use '
                                   '[default = 1]')
    optimisation.add_argument('--device-id',
                              type=int,
                              default=0,
                              help='Device ID of GPU to use '
                                   '[default = 0]')

    other = parser.add_argument_group('Other')
    other.add_argument('--version', action='version',
                       version='%(prog)s '+__version__)

    return parser.parse_args()


def main():
    args = get_options()

    kmers = sorted([int(i) for i in args.k.split(",")])
    for k in kmers:
        if k % 2 == 0:
            raise RuntimeError("Even k-mer lengths not allowed")

    if args.assemble:
        pass

        sparrowhawk.assembleReads(args.fwd,
                                    args.rev,
                                    kmers,
                                    output,
                                    args.cpus,
                                    args.device_id)

    elif args.dbg:
        pass

        if len(kmers) > 1:
            sys.stderr.write("Using smallest k-mer only\n")
        k = kmers[0]

        sparrowhawk.makeDBG(args.fwd,
                             args.rev,
                             k,
                             output,
                             args.cpus,
                             args.device_id)

    elif args.info:
        sparrowhawk_py = __version__
        sparrowhawk_cu = sparrowhawk.version
        cuda_version = sparrowhawk.cudaVersion
        print("sparrowhawk python module version: " + sparrowhawk_py)
        print("sparrowhawk CUDA module version: " + sparrowhawk_cu + \
              " (CUDA v" + cuda_version + ")")

        print("GPUs available:")
        print(sparrowhawk.deviceInfo())

    sys.exit(0)

if __name__ == "__main__":
    main()