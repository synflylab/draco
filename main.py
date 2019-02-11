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
parser.add_argument('sequence', help='TALE binding sequence (DNA)')
parser.add_argument('--upstream', help='Protein sequence to include upstream of the repeats', default='', type=str)
parser.add_argument('--downstream', help='Protein sequence to include downstream of the repeats', default='', type=str)
parser.add_argument('--avoid', action='append', help='Avoid this fragment in the optimized sequence', type=str)
parser.add_argument('--check-repeats', help='Check the optimized sequence for repeats', default=True, type=bool)
parser.add_argument('--repeat-len', default=20, help='Max allowed length of a repeat', type=int)
parser.add_argument('--check-inv-repeats', help='Check the optimized sequence for inverse repeats', default=True, type=bool)
parser.add_argument('--inv-repeat-len', default=12, help='Max allowed length of an inverted repeat', type=int)
parser.add_argument('--check-stretch', help='Check for base stretches', default=True, type=bool)
parser.add_argument('--stretch-len', default=8, help='Max allowed length of a base stretch', type=int)
parser.add_argument('--check-gc', help='Check GC content', default=True, type=bool)
parser.add_argument('--max-gc', default=65, help='Max allowed GC content', type=int)
parser.add_argument('--gc-window', default=200, help='GC calculation window', type=int)
parser.add_argument('--repeat-attempts', default=100, help='Max attempts to optimize a repeat', type=int)
parser.add_argument('--seq-attempts', default=1000, help='Max attempts to optimize the sequence', type=int)
parser.add_argument('--handicap', default=1, help='Max codon handicap for the optimized sequence', type=int)
parser.add_argument('--progress', help='Show progress indicator', action='store_true', type=bool)
parser.add_argument('--log', help='Logging level')
parser.add_argument('--debug', help='Enable debug logging')
args = parser.parse_args()

upstream = 'DTGQLVKIAKRGGVTAMEAVHASRNALTGAPLN'
downstream = 'SIVAQLSRPDPALAALTNDHLVALACLGGRPAM'

if args.debug:
    args.log = 'DEBUG'

if args.log:
    logging.basicConfig(level=args.log.upper())
    logging.getLogger('BiopythonWarning').setLevel(logging.INFO)

target = Seq(args.sequence, alphabet=IUPAC.ambiguous_dna)
upstream = Seq(args.upstream, alphabet=IUPAC.protein)
downstream = Seq(args.downstream, alphabet=IUPAC.protein)
tale = TALE(target, upstream=upstream, downstream=downstream)

codons = CodonUsage(CodonUsage.Dmel)
draco = Draco(tale, codons, args.avoid,
              args.check_repeats, args.repeat_len, args.repeat_len,
              args.check_inv_repeats, args.inv_repeat_len, args.inv_repeat_len,
              args.check_stretch, args.stretch_len,
              args.check_gc, args.max_gc, args.gc_window,
              args.repeat_attempts, args.seq_attempts,
              args.handicap, args.progress)
dna = draco.random().back_transcribe()
if dna:
    print('\n' + str(dna))
