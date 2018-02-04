#!/usr/bin/env python

# Author: Aretas Gaspariunas

###PLEASE READ!###
## Detects ORFs in a DNA sequence provided as input FASTA file
## Assumes that the starting codon is ATG, accounts all 6 reading frames
## By default ignores ORFs containing ambigious calls indicated by letter 'N'
## Does not ignore nested ORFs
## Please use the help flag for further information and options
###PLEASE READ!###

import re
from Bio.Seq import Seq  # Biopython module; install with pip via CLI

class DNAsequence:
    def __init__(self, sequence, full_tag):
        self.sequence = sequence.upper()
        self.full_tag = full_tag
        self.tag = full_tag.split(' ')[0]

    # a method to find a total number of a specified codon
    find_codon_len = lambda self, codon: len(re.findall(codon, self.sequence))

    # a method to find GC precentage in the sequence
    def count_GC(self):
        g_count = self.sequence.count('G')
        c_count = self.sequence.count('C')
        gc_count = float((g_count + c_count)) / float(len(self.sequence)) * 100

        return gc_count

    # a method to find a reverse compliment strand
    def rev_compliment(self):
        translate_dict = {ord('A'):'T', ord('T'):'A', ord('C'):'G', ord('G'):'C'}
        return self.sequence.translate(translate_dict)[::-1]

    # a method to translate a DNA sequence to an amino acid sequence
    def translate(self, table):
        return Seq(self.sequence).translate(table)

    # a method to find ORFs in a sequence (3 reading frames)
    def find_ORF (self, seq, min_nt_len, strand, start_codon, nan):
        orf_dict = {}
        frame = 0
        while frame < 3:
            for i in range(frame, len(seq), 3):
                if seq[i:i+3] == '{0}'.format(start_codon):
                    for i2 in range(i+3, len(seq), 3):
                        if seq[i2:i2+3] in ['TAG', 'TAA', 'TGA']:
                            if len(seq[i:i2+3]) >= min_nt_len:
                                if nan is False and re.search('N', seq[i:i2+3]):
                                    print ("Ignoring ORF starting at {0}; contains an ambigious NT call: 'N'".format(i+1))
                                    pass
                                elif i2+3 not in orf_dict.keys():
                                    orf_dict[i2+3] = CodingDNA(seq[i:i2+3], self.full_tag, strand, frame + 1, i + 1) # self.sequence[i:i2+3]
                                elif i2+3 in orf_dict.keys():
                                    if len(orf_dict[i2+3].sequence) > len(seq[i:i2+3]):
                                        pass
                            break
            frame += 1

        return orf_dict

class CodingDNA(DNAsequence):
    def __init__(self, sequence, full_tag, strand, frame, start):
        super().__init__(sequence, full_tag)
        self.strand = strand
        self.frame = frame
        self.start = start

if __name__ == '__main__':
    pass
