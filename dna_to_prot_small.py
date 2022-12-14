#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: asansal
Using partial code from:
Prestevez - dna2proteins: https://github.com/prestevez/dna2proteins
Copyright (c) 2015 Patricio Rodrigo Estévez Soto
brentp - gistfile1.py: https://gist.github.com/brentp/477969#file-gistfile1-py
"""


NUCLEOTIDE_BASE = {
    "DNA": ["A", "T", "C", "G"],
    "RNA": ["A", "U", "C", "G"]
}

DNA_Codons = {
    # 'M' - START, '_' - STOP
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "_", "TAG": "_", "TGA": "_"
}

RNA_Codons = {
    # 'M' - START, '_' - STOP
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UGU": "C", "UGC": "C",
    "GAU": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "UUU": "F", "UUC": "F",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAU": "H", "CAC": "H",
    "AUA": "I", "AUU": "I", "AUC": "I",
    "AAA": "K", "AAG": "K",
    "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUG": "M",
    "AAU": "N", "AAC": "N",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UGG": "W",
    "UAU": "Y", "UAC": "Y",
    "UAA": "_", "UAG": "_", "UGA": "_"
}

def read_FASTA(filePath):
    with open(filePath, 'r') as f:
        FASTAFile = [l.strip() for l in f.readlines()]

    FASTADict = {}
    FASTALabel = ""

    for line in FASTAFile:
        if '>' in line:
            FASTALabel = line
            FASTADict[FASTALabel] = ""
        else:
            FASTADict[FASTALabel] += line

    return FASTADict

import pandas
import sys
from Bio import SeqIO
import csv

class bio_seq:
    """DNA sequnece class. Defalt value: ATCG, DNA, No label"""

    def __init__(self, seq="ATCG", seq_type="DNA", label='No Label'):
        """Sequence initialization, validation."""
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.__validate()
        assert self.is_valid, f"Provided data does not seem to be a correct {self.seq_type} sequence"

    # DNA Toolkit functions:
    
    def __validate(self):
        """Check the sequence to make sure it is a valid DNA string"""
        return set(NUCLEOTIDE_BASE[self.seq_type]).issuperset(self.seq)

    def transcription(self):
        """DNA -> RNA Transcription. Replacing Thymine with Uracil"""
        if self.seq_type == "DNA":
            return self.seq.replace("T", "U")
        return "Not a DNA sequence"

    def reverse_complement(self):
        """
        Swapping adenine with thymine and guanine with cytosine.
        Reversing newly generated string
        """
        if self.seq_type == "DNA":
            mapping = str.maketrans('ATCG', 'TAGC')
        else:
            mapping = str.maketrans('AUCG', 'UAGC')
        return self.seq.translate(mapping)[::-1]

    def translate_seq(self, init_pos=0):
        """Translates a DNA sequence into an aminoacid sequence"""
        if self.seq_type == "DNA":
            return [DNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]
        elif self.seq_type == "RNA":
            return [RNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]
    
    def gen_reading_frames_fw(self):
        """Generate the six reading frames of a DNA sequence, including reverse complement"""
        frames = []
        frames.append(self.translate_seq(0))
        frames.append(self.translate_seq(1))
        frames.append(self.translate_seq(2))
        return frames
    
    def gen_reading_frames_rv(self):
        frames = []
        tmp_seq = bio_seq(self.reverse_complement(), self.seq_type)
        frames.append(tmp_seq.translate_seq(0))
        frames.append(tmp_seq.translate_seq(1))
        frames.append(tmp_seq.translate_seq(2))
        del tmp_seq
        return frames
    
    def proteins_from_rf(self, aa_seq, init, seqname, sense = 'fw'):
        """Compute all possible proteins in an aminoacid seq and return a list of possible proteins.
        sense = ('fw', 'rv')
        """
        current_prot, coordinates_prot, proteins, prot_coord = [],[],[],[]
        
        temp_init=-3+init
        endposition=-1+init
        initposition=0
        
        for aa in aa_seq:
            endposition += 3
            temp_init += 3
            
            n=len(current_prot)
            
            if aa == "_":
                # STOP accumulating amino acids if _ - STOP was found
                if current_prot:
                    
                    if n > 10:
                        #Length between 11 and 100
                        if n < 101:
                            proteins.append(current_prot)
                            prot=''.join(current_prot)

                            if sense == 'fw':
                                string=f'{symbol}{seqname}{sep1}{initposition}{sep2}{endposition}'
                                coordinates_prot= tuple([initposition, endposition, prot, string])
                                prot_coord.append(coordinates_prot)

                            elif sense == 'rv':
                                coordinates_prot= tuple([initposition, endposition, prot])
                                prot_coord.append(coordinates_prot)

                            else:
                                print('Only options are "fw" or "rv"')

                            current_prot = []
                            coordinates_prot = []
                        else:
                            current_prot = []
                    else:
                        current_prot = []
                   
                        
            else:
                # START accumulating amino acids if M - START was found
                if aa == "M":
                    #Only if this M is the 1st, not the following ones
                    if n == 0:
                        initposition=temp_init
                        current_prot=['M']
                        
                if n > 0:
                    #If protein already iniciated, add aa
                    current_prot += aa

        return prot_coord

'''
The usage will be:
    python thiscode.py [1]fastafile [2]wanted [3]separation
'''

#If we don't write 'all' we'll get the names from wanted.txt
if sys.argv[2] != 'all':
    wanted = [line.strip() for line in open(sys.argv[2])]
#If it's 'all' written, we'll get ALL the sequences from the fasta file
else:
    wanted=list(SeqIO.index(sys.argv[1], "fasta"))

#Separation has 2 options: 'joined' or 'separated'
separation = sys.argv[3]

DNAsequenceFASTA = read_FASTA(sys.argv[1])

#Format of the initial fasta file
symbol = '>'

#Format of the 'string' file
sep1, sep2, rv = ':', '..', ':r'

#format of the output file
format_val='.prot.fasta'

joined_data=pandas.DataFrame()

for seq in wanted:
    
    DNAsequence = bio_seq(DNAsequenceFASTA[f'{symbol}{seq}'], seq_type="DNA")
    seqname = seq
    
    #Total nucleotides number for calculating coordinates
    DNAlength=len(DNAsequence.seq) 
    
    rfs_fw = DNAsequence.gen_reading_frames_fw()
    rfs_rv = DNAsequence.gen_reading_frames_rv()
    
    res = []
    init=0
    
    #Loop for forward sequences
    for rf in rfs_fw:
        init+=1
        prots = DNAsequence.proteins_from_rf(rf, init, seqname, 'fw')
        for p in prots:
            res.append(p)
    
    #Loop for reverse sequences
    init=-1
    for rf in rfs_rv:
        init+=1
        prots = DNAsequence.proteins_from_rf(rf, init, seqname, 'rv')
        for p in prots:
            #Obtaining coordinates from reverse sequences
            reverse_p1=DNAlength-p[0]
            reverse_p2=DNAlength-p[1]
            
            #Name of the fasta sequence
            string=f'{symbol}{seqname}{sep1}{reverse_p2}{sep2}{reverse_p1}{rv}'
            
            reverse_p=tuple([reverse_p2,reverse_p1, p[2], string])
            res.append(reverse_p)
    
    df = pandas.DataFrame(res)
    df = df.sort_values(0)
    subset_prot = df[2]
    subset_name = df[3]
    
    fasta_prot = pandas.DataFrame()
    
    #Extracting data from the df
    for row in range(len(df)):
        name_prot = pandas.DataFrame([subset_name[row]])
        sequence_prot = pandas.DataFrame([subset_prot[row]])
        new_prot = pandas.concat([name_prot,sequence_prot], ignore_index=True, sort=False, axis=0)
        fasta_prot = pandas.concat([fasta_prot, new_prot], ignore_index=True, sort=False, axis=0)
    
    #If joined, we save it all in the same variable
    if separation == 'joined':
        joined_data = pandas.concat([joined_data, fasta_prot], ignore_index=True, sort=False, axis=0)
    
    #If separated, we extract a new one each time in the loop
    elif separation == 'separated':
        fasta_prot.to_csv(f'{seq}{format_val}', sep = '\n', header=False, index=False)
    
    else:
        print("Only options: 'joined' or 'separated'")
        
#Finally we extract the fasta file for the joined sequences
seq = 'joined'
if separation == 'joined':
    joined_data.to_csv(f'{seq}{format_val}', sep = '\n', line_terminator='\n', header=False, index=False, 
                       quoting=csv.QUOTE_NONE, quotechar="",  escapechar="\\")




