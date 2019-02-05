#!/usr/bin/env python3

import argparse
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from tale import TALE
from draco import Draco
from codon import CodonUsage

parser = argparse.ArgumentParser(description='Direct Repeat Aware Codon Optimizer')
parser.add_argument('sequence', help='TALE binding sequence')
args = parser.parse_args()

print(args.sequence)
target = Seq(args.sequence, alphabet=IUPAC.ambiguous_dna)
tale = TALE(args.sequence)
print(tale)

codons = CodonUsage(CodonUsage.dmel)
print(codons)

print(codons['CUC'].fraction)
print(codons['L'])
