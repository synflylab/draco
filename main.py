#!/usr/bin/env python3

import argparse
import logging
import warnings
from Bio import BiopythonWarning
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from tale import TALE
from draco import Draco
from codon import CodonUsage

warnings.simplefilter('ignore', BiopythonWarning)

parser = argparse.ArgumentParser(description='Direct Repeat Aware Codon Optimizer')
parser.add_argument('sequence', help='TALE binding sequence')
parser.add_argument('--check-repeats', help='Check the optimized sequence for repeats', default=True)
parser.add_argument('--repeat-len', default=20, help='Max allowed length of a repeat')
parser.add_argument('--check-inv-repeats', help='Check the optimized sequence for inverse repeats', default=True)
parser.add_argument('--inv-repeat-len', default=12, help='Max allowed length of an inverted repeat')
parser.add_argument('--repeat-attempts', default=100, help='Max attempts to optimize a repeat')
parser.add_argument('--seq-attempts', default=1000, help='Max attempts to optimize the sequence')
parser.add_argument('--handicap', default=1, help='Max codon handicap for the optimized sequence')
parser.add_argument('--progress', help='Show progress indicator', action='store_true')
parser.add_argument('--log', help='Logging level')
parser.add_argument('--debug', help='Enable debug logging', action='store_true')
args = parser.parse_args()

if args.debug:
    args.log = 'DEBUG'

if args.log:
    logging.basicConfig(level=args.log.upper())
    logging.getLogger('BiopythonWarning').setLevel(logging.INFO)

target = Seq(args.sequence, alphabet=IUPAC.ambiguous_dna)
tale = TALE(args.sequence)

codons = CodonUsage(CodonUsage.Dmel)
draco = Draco(tale, codons, args.check_repeats, args.repeat_len, args.repeat_len, args.check_inv_repeats,
              args.inv_repeat_len, args.inv_repeat_len, args.repeat_attempts, args.seq_attempts, args.handicap)
dna = draco.random().back_transcribe()
if dna:
    print('\n' + str(dna))
