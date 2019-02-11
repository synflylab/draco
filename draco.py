import random
import re
import logging
from Bio.Alphabet import IUPAC
from Bio.Data import IUPACData
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from progress.bar import Bar

class Draco:
    """
    DRACO - Direct Repeat Aware Codon Optimizer
    This class computes the codon optimized repeat sequence
    """

    min_repeat_len = 20
    max_repeat_len = 20
    min_invrep_len = 12
    max_invrep_len = 12
    max_stretch = 8
    max_gc = 65
    gc_window = 200
    max_repeat_attempts = 100
    max_sequence_attempts = 1000
    max_handicap = 1

    def __init__(self, sequence, codons, avoid,
                 check_repeats=True, min_repeat_len=None, max_repeat_len=None,
                 check_invreps=True, min_invrep_len=None, max_invrep_len=None,
                 check_stretch=True, max_stretch=None,
                 check_gc=True, max_gc=None, gc_window=None,
                 max_repeat_attempts=None, max_sequence_attempts=None,
                 max_handicap=None, progress=True):

        self.sequence = sequence
        self.codons = codons
        self.avoid = avoid
        self.check_repeats = check_repeats
        self.check_invreps = check_invreps
        self.check_stretch = check_stretch
        self.check_gc = check_gc
        self.min_repeat_len = self._default(min_repeat_len, self.min_repeat_len)
        self.max_repeat_len = self._default(max_repeat_len, self.max_repeat_len)
        self.min_invrep_len = self._default(min_invrep_len, self.min_invrep_len)
        self.max_invrep_len = self._default(max_invrep_len, self.max_invrep_len)
        self.max_stretch = self._default(max_stretch, self.max_stretch)
        self.max_gc = self._default(max_gc, self.max_gc)
        self.gc_window = self._default(gc_window, self.gc_window)
        self.max_repeat_attempts = self._default(max_repeat_attempts, self.max_repeat_attempts)
        self.max_sequence_attempts = self._default(max_sequence_attempts, self.max_sequence_attempts)
        self.max_handicap = self._default(max_handicap, self.max_handicap)
        self.progress = progress
        self.bar = None
        self.residues = {}
        self.variants = {}
        self.max_ocr = {}
        self.repeat_ocr = {}
        self.thresholds = {}

    @staticmethod
    def _default(value, default):
        return default if value is None else value

    def _min_max_len(self, min_length, max_length, inverse):
        if inverse:
            min_length = self._default(min_length, self.min_invrep_len)
            max_length = self._default(max_length, self.max_invrep_len)
        else:
            min_length = self._default(min_length, self.min_repeat_len)
            max_length = self._default(max_length, self.max_repeat_len)

        return min_length, max_length

    def _analyse_codons(self, handicap=0):
        """
        Analyze codon distribution for the repeat sequence

        :param handicap:
        :return:
        """
        self.residues = {}
        self.variants = {}
        self.max_ocr = {}
        self.repeat_ocr = {}
        self.thresholds = {}
        sequence = self.sequence.sequence()
        n_repeats = len(self.sequence.repeats)
        for r in sequence:
            if r not in IUPACData.protein_letters:
                continue
            if r not in self.residues.keys():
                self.residues[r] = 0
                self.variants[r] = len(self.codons[r])
            self.residues[r] = self.residues[r] + 1

        for r in self.residues.keys():
            self.max_ocr[r] = {}
            self.repeat_ocr[r] = {}
            for codon in self.codons[r]:
                self.max_ocr[r][codon.rna] = round(float(codon.fraction) * self.residues[r]) + handicap
                self.repeat_ocr[r][codon.rna] = self.max_ocr[r][codon.rna] / n_repeats

    def random(self):
        """
        Use guided random to create codon-optimized sequence

        :return: string RNA sequence
        """
        count = 1
        handicap = 0
        if self.progress:
            self.bar = Bar()
        while True:
            if self.progress:
                self.bar.index = 0
                self.bar.message = 'Computing sequence [' + str(count) + '/' + str(self.max_sequence_attempts) + ']'
            logging.debug("Computing sequence: " + str(count) + " attempt...")
            self._analyse_codons(handicap)
            sequence = self._compute_sequence()
            if sequence is not False:
                logging.debug("Success!")
                break
            count += 1
            if count % (self.max_sequence_attempts / (self.max_handicap + 1)) == 0:
                handicap += 1
                logging.info("Increasing handicap to", handicap)
            if count == self.max_sequence_attempts:
                logging.error("Failed to find a suitable sequence! Try increasing the number of attempts,"
                              " the handicap or max repeat length")
                return ""
        return sequence

    def _get_repeat_words(self, repeat_a, repeat_b, min_length=None, max_length=None, inverse=False):
        """
        Retrieve words from a pair of consecutive repeats

        :param repeat_a:
        :param repeat_b:
        :param min_length:
        :param max_length:
        :param inverse:
        :return:
        """

        min_length, max_length = self._min_max_len(min_length, max_length, inverse)
        sequence = repeat_a + repeat_b
        end = len(sequence)
        for length in range(max_length, min_length - 1, -1):
            start = 0 if len(repeat_a) - length < 0 else len(repeat_a) - length
            for pos in range(start, end - min_length):
                word = sequence[pos:(pos + length)]
                if inverse:
                    yield word.reverse_complement()
                else:
                    yield word

    def _get_words(self, sequence, min_length=None, max_length=None, inverse=False):
        """
        Retrieve words from a generic sequence

        :param sequence:
        :param min_length:
        :param max_length:
        :param inverse:
        :return:
        """

        min_length, max_length = self._min_max_len(min_length, max_length, inverse)
        end = len(sequence)
        for start in range(0, end - min_length):
            for length in range(max_length, min_length - 1, -1):
                word = sequence[start:(start + length)]
                if inverse:
                    yield word.reverse_complement()
                else:
                    yield word

    def _find_repeats(self, sequence, min_length=None, max_length=None, words=None, inverse=False):
        """
        Find direct repeats in a sequence

        :param sequence:
        :param min_length:
        :param max_length:
        :param words:
        :return:
        """

        min_length, max_length = self._min_max_len(min_length, max_length, inverse)
        repeats = {}
        dups = {}
        if words is None:
            words = self._get_words(sequence, min_length, max_length, inverse)
        for word in words:
            self._find_all(sequence, word, start=0, repeats=repeats, duplicates=dups)
        return repeats

    @staticmethod
    def _find_all(sequence, word, start=0, end=0, repeats=None, duplicates=None):
        """
        Find all occurrences of a word in sequence

        :param sequence:
        :param word:
        :param start:
        :param end:
        :param repeats:
        :param duplicates:
        :return:
        """
        start = start + 1
        end = len(sequence) if end == 0 else end
        repeats = {} if repeats is None else repeats
        duplicates = {} if duplicates is None else duplicates
        length = len(word)
        matches = []
        count = 0
        while True:
            start = sequence.find(word, start, end)
            if start == -1:
                break
            matches.append(start)
            start += length
        if len(matches) > 1:
            curr = tuple(matches)
            prev = tuple([x - 1 for x in matches])
            if curr not in duplicates.keys() and prev not in duplicates.keys():
                repeats[curr] = word
                count += 1
            duplicates[curr] = length
        return count

    def _check_fragment(self, rna, sequences, index):
        """
        Check if sequence has repeats
        """
        repeats = []
        inv_repeats = []
        stretches = []
        gc = 0

        seq = sequences[index - 1] + rna

        if self.check_repeats:
            words = self._get_repeat_words(sequences[index - 1], rna)
            repeats = self._find_repeats(sum(sequences, Seq('', IUPAC.ambiguous_rna)) + rna, words=words)
        if self.check_invreps:
            inv_words = self._get_repeat_words(sequences[index - 1], rna, inverse=True)
            inv_repeats = self._find_repeats(sum(sequences, Seq('', IUPAC.ambiguous_rna)) + rna,
                                             words=inv_words, inverse=True)
        if self.check_stretch:
            stretches = re.findall(r'((\w)\2{' + str(self.max_stretch - 1) + ',})', str(seq))
        if self.check_gc:
            end = len(seq) - self.gc_window if len(seq) > self.gc_window else 1
            for start in range(0, end):
                cur_gc = GC(seq[start:start+self.gc_window])
                gc = cur_gc if cur_gc > gc else gc
                if gc > self.max_gc:
                    break

        forbidden = False
        if len(self.avoid > 0):
            for word in self.avoid:
                forbidden = (seq.find(word) != -1)
                if forbidden:
                    break

        return len(repeats) == 0 and len(inv_repeats) == 0 and len(stretches) == 0 \
            and gc < self.max_gc and (not forbidden)

    def _compute_sequence(self):
        """
        Compute codon-optimized repeat sequence

        :return:
        """

        fragments = []
        if self.sequence.upstream:
            fragments.append(self.sequence.upstream)
        fragments += [r.sequence() for r in self.sequence.repeats]
        if self.sequence.downstream:
            fragments.append(self.sequence.downstream)

        sequences = []
        if self.progress:
            self.bar.max = len(fragments)
        for index, fragment in enumerate(fragments):
            take = 0
            while True:
                max_ocr = self.max_ocr.copy()
                residues = self.residues.copy()
                rna = self._compute_fragment(fragment)
                if index > 0:
                    if self._check_fragment(rna, sequences, index):
                        break
                    if take > self.max_repeat_attempts:
                        return False
                else:
                    break
                self.max_ocr = max_ocr
                self.residues = residues
                take += 1
            sequences.append(rna)
            if self.progress:
                self.bar.next()

        return sum(sequences, Seq('', IUPAC.ambiguous_rna))

    def _compute_fragment(self, sequence):
        """
        Codon-optimize a single repeat

        :param sequence:
        :return:
        """

        rna = Seq('', alphabet=IUPAC.ambiguous_rna)
        k = 0
        for r in sequence:
            if r not in IUPACData.protein_letters:
                continue
            self._recalculate_thresholds(r)
            toss = random.randint(0, self.residues[r])
            for k, t in enumerate(self.thresholds[r]):
                if toss <= t:
                    break
            codon = self.codons[r][k]
            rna = rna + codon.rna
            self.max_ocr[r][codon.rna] = self.max_ocr[r][codon.rna] - 1
            self.residues[r] = self.residues[r] - 1
        return rna

    def _recalculate_thresholds(self, r):
        """
        Recalculate codon-usage thresholds

        :param r:
        :return:
        """

        total = 0
        self.thresholds[r] = []
        for codon in self.max_ocr[r]:
            total = total + self.max_ocr[r][codon]
            self.thresholds[r].append(total)
