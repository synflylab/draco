from Bio.Alphabet import IUPAC
from Bio.Data import IUPACData
from Bio.Seq import Seq

class CodonUsage:
    Dmel = """
UUU F 0.38 13.2 (289916)  UCU S 0.08  7.0 (154186)  UAU Y 0.37 10.8 (236811)  UGU C 0.29  5.4 (118088)
UUC F 0.62 21.8 (479372)  UCC S 0.24 19.6 (429341)  UAC Y 0.63 18.4 (403675)  UGC C 0.71 13.2 (288853)
UUA L 0.05  4.5 (097715)  UCA S 0.09  7.8 (171695)  UAA * 0.41  0.8 (017807)  UGA * 0.25  0.5 (010767)
UUG L 0.18 16.1 (353621)  UCG S 0.20 16.6 (365159)  UAG * 0.33  0.7 (014362)  UGG W 1.00  9.9 (217518)

CUU L 0.10  9.0 (196787)  CCU P 0.13  6.9 (151856)  CAU H 0.40 10.8 (236061)  CGU R 0.16  8.8 (192276)
CUC L 0.15 13.8 (303153)  CCC P 0.33 18.1 (396168)  CAC H 0.60 16.2 (354699)  CGC R 0.33 18.0 (395106)
CUA L 0.09  8.2 (180360)  CCA P 0.25 13.5 (297071)  CAA Q 0.30 15.6 (342415)  CGA R 0.15  8.4 (185119)
CUG L 0.43 38.2 (839127)  CCG P 0.29 15.8 (347206)  CAG Q 0.70 36.1 (792657)  CGG R 0.15  8.2 (180473)

AUU I 0.34 16.6 (363497)  ACU T 0.17  9.5 (208889)  AAU N 0.44 21.0 (460669)  AGU S 0.14 11.5 (252555)
AUC I 0.47 22.9 (502821)  ACC T 0.38 21.3 (467509)  AAC N 0.56 26.2 (575297)  AGC S 0.25 20.4 (447808)
AUA I 0.19  9.5 (208315)  ACA T 0.20 11.0 (241893)  AAA K 0.30 17.0 (372524)  AGA R 0.09  5.1 (112784)
AUG M 1.00 23.6 (518200)  ACG T 0.26 14.4 (315479)  AAG K 0.70 39.5 (866960)  AGG R 0.11  6.3 (137902)

GUU V 0.19 11.0 (240735)  GCU A 0.19 14.4 (315879)  GAU D 0.53 27.6 (604730)  GGU G 0.21 13.3 (291161)
GUC V 0.24 13.9 (304893)  GCC A 0.45 33.6 (736394)  GAC D 0.47 24.6 (540386)  GGC G 0.43 26.7 (587016)
GUA V 0.11  6.4 (139476)  GCA A 0.17 12.8 (280181)  GAA E 0.33 21.1 (462468)  GGA G 0.29 18.0 (395377)
GUG V 0.47 27.8 (609794)  GCG A 0.19 14.0 (307977)  GAG E 0.67 42.5 (933622)  GGG G 0.07  4.7 (102708)
"""

    def __init__(self, string):
        self.codons = {}
        self.residues = {}
        self._parse(string)

    def __getitem__(self, key):
        if key in self.codons.keys():
            return self.codons.get(key)
        elif key in self.residues.keys():
            return self.residues.get(key)
        else:
            raise KeyError

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        string = ""
        separator = ""
        count = 0
        for key, codon in self.codons.items():
            string = string + separator + codon.rna + '(' + codon.protein + ')' +\
                     ':' + "{: >6.2f}".format(float(codon.fraction))
            count = count + 1
            if count % 4 == 0:
                separator = "\n"
            else:
                separator = "  "
        return string

    def _parse(self, string):
        for line in string.splitlines():
            definitions = line.split(')  ')
            for definition in definitions:
                if definition.strip():
                    codon = Codon(definition)
                    self.codons[codon.rna] = codon
                    if codon.protein not in self.residues.keys():
                        self.residues[codon.protein] = []
                    self.residues[codon.protein].append(codon)


class Codon:
    def __init__(self, string, protein=None, fraction=None, frequency=None):
        if protein is None:
            parsed = string.split()
            self.rna = Seq(parsed[0], alphabet=IUPAC.ambiguous_rna)
            self.protein = Seq(parsed[1], alphabet=IUPAC.protein)
            self.fraction = parsed[2]
            self.frequency = parsed[3]
        else:
            self.rna = Seq(string, alphabet=IUPAC.ambiguous_rna)
            self.protein = Seq(protein, alphabet=IUPAC.protein)
            self.fraction = fraction
            self.frequency = frequency

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return self.rna + '(' + self.protein + '): ' + "{:.2f}".format(float(self.fraction))
