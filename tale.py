from Bio.Alphabet import IUPAC
from Bio.Data import IUPACData
from Bio.Seq import Seq

class TALrepeat:

    A = 'LTPDQVVAIAS'
    B = 'GGKQALETVQRLLPVLCQDHG'
    C = 'GGKQALE'

    def __init__(self, nucleotide, last=False):
        self.last = last
        if nucleotide.upper() == 'G' or nucleotide.upper() == 'R':
            self.rvd = 'NN'
        elif nucleotide.upper() == 'A':
            self.rvd = 'NI'
        elif nucleotide.upper() == 'T':
            self.rvd = 'NG'
        elif nucleotide.upper() == 'C':
            self.rvd = 'HD'
        elif nucleotide in IUPACData.ambiguous_dna_letters:
            self.rvd = 'NS'
        else:
            self.rvd = None

    def __str__(self):
        return self.sequence()

    def sequence(self):
        if self.rvd is None:
            return None
        else:
            if self.last:
                return Seq(self.A + self.rvd + self.C, alphabet=IUPAC.protein)
            else:
                return Seq(self.A + self.rvd + self.B, alphabet=IUPAC.protein)

    def template(self):
        return Seq(self.A + 'XX' + self.B, alphabet=IUPAC.protein)


class TALE:

    def __init__(self, target, upstream=None, downstream=None):
        self.target = target
        self.repeats = []
        self.upstream = upstream
        self.downstream = downstream
        index = 0
        for n in target:
            index += 1
            if index == len(target):
                self.repeats.append(TALrepeat(n, True))
            else:
                self.repeats.append(TALrepeat(n))

    def __str__(self):
        return self.protein()

    def protein(self):
        sequence = self.upstream
        for r in self.repeats:
            sequence += r.sequence()
        sequence += self.downstream

        return sequence

    def sequence(self):
        return self.protein()

    def repeat(self):
        return self.repeats[0].template()

    def n_repeats(self):
        return len(self.repeats)
