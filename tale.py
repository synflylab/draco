from Bio.Data import IUPACData

class TALrepeat:

    A = 'LTPDQVVAIAS'
    B = 'GGKQALETVQRLLPVLCQDHG'

    def __init__(self, nucleotide):
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
            return self.A + self.rvd + self.B


class TALE:

    def __init__(self, target):
        self.target = target
        self.repeats = []
        for n in target:
            self.repeats.append(TALrepeat(n))

    def __str__(self):
        return self.sequence()

    def sequence(self):
        sequence = ''
        for r in self.repeats:
            sequence = sequence + r.sequence() + '\n'

        return sequence
