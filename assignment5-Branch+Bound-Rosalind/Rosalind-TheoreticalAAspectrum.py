#!/usr/bin/env python3
########################################################################
# File:Rosalind-TheoreticalAAspectrum.py
#  executable: Rosalind-TheoreticalAAspectrum.py <input >output
# Purpose: Find the Theoretical spectrum of a given amino acid string.
#   stdout: Theoretical spectrum of amino acid sequence
#
# Author:       Akshar Lohith
# Group(sounding boards): Lauren Sanders, Roger Volden
# History:      AL 10/29/2016 Created
#
#Program flow method
########################################################################
import sys
class CycloPeptide:
    """docstring for CycloPeptide."""
    # Single aminoAcid to interger mass
    aa2IntMass = {'G':57,'A':71,'S':87,'P':97,'V':99,'T':101,'C':103,'I':113,'L':113,
              'N':114,'D':115,'K':128,'Q':128,'E':129,'M':131,'H':137,'F':147,
              'R':156,'Y':163,'W':186}
    # As written, these are accessed as class attributes, for example:
    # ProteinParam.aa2mw['A'] or ProteinParam.mwH2O
    # Taken from BME160 ProteinParams class.
    RNAcodonTable = {
    #                        Second Base
    #        U             C             A             G
    # U
        'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys',     # UxU
        'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys',     # UxC
        'UUA': 'Leu', 'UCA': 'Ser', 'UAA': '---', 'UGA': '---',     # UxA
        'UUG': 'Leu', 'UCG': 'Ser', 'UAG': '---', 'UGG': 'Trp',     # UxG
    # C
        'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg',     # CxU
        'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg',     # CxC
        'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg',     # CxA
        'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',     # CxG
    # A
        'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser',     # AxU
        'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser',     # AxC
        'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg',     # AxA
        'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',     # AxG
    # G
        'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly',     # GxU
        'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly',     # GxC
        'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly',     # GxA
        'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'      # GxG
                      }

    DNAcodonTable = {
        #table made based off of Model's RNA condon table
        'TTT':'Phe', 'TCT':'Ser', 'TAT':'Tyr', 'TGT':'Cys',
        'TTC':'Phe', 'TCC':'Ser', 'TAC':'Tyr', 'TGC':'Cys',
        'TTA':'Leu', 'TCA':'Ser', 'TAA':'---', 'TGA':'---',
        'TTG':'Leu', 'TCG':'Ser', 'TAG':'---', 'TGG':'Trp',
        'CTT':'Leu', 'CCT':'Pro', 'CAT':'His', 'CGT':'Arg',
        'CTC':'Leu', 'CCC':'Pro', 'CAC':'His', 'CGC':'Arg',
        'CTA':'Leu', 'CCA':'Pro', 'CAA':'Gln', 'CGA':'Arg',
        'CTG':'Leu', 'CCG':'Pro', 'CAG':'Gln', 'CGG':'Arg',
        'ATT':'Ile', 'ACT':'Thr', 'AAT':'Asn', 'AGT':'Ser',
        'ATC':'Ile', 'ACC':'Thr', 'AAC':'Asn', 'AGC':'Ser',
        'ATA':'Ile', 'ACA':'Thr', 'AAA':'Lys', 'AGA':'Arg',
        'ATG':'Met', 'ACG':'Thr', 'AAG':'Lys', 'AGG':'Arg',
        'GTT':'Val', 'GCT':'Ala', 'GAT':'Asp', 'GGT':'Gly',
        'GTC':'Val', 'GCC':'Ala', 'GAC':'Asp', 'GGC':'Gly',
        'GTA':'Val', 'GCA':'Ala', 'GAA':'Glu', 'GGA':'Gly',
        'GTG':'Val', 'GCG':'Ala', 'GAG':'Glu', 'GGG':'Gly'
                      }
    aasingle2triple = {
        'G':'Gly', 'A':'Ala', 'V':'Val', 'L':'Leu', 'I':'Ile',
        'M':'Met', 'F':'Phe', 'W':'Trp', 'P':'Pro', 'S':'Ser',
        'T':'Thr', 'C':'Cys', 'Y':'Tyr', 'N':'Asn', 'Q':'Gln',
        'D':'Asp', 'E':'Glu', 'K':'Lys', 'R':'Arg', 'H':'His'
                      }
    aatriple2single = {
        'GLY':'G', 'ALA':'A', 'VAL':'V', 'LEU':'L', 'ILE':'I',
        'MET':'M', 'PHE':'F', 'TRP':'W', 'PRO':'P', 'SER':'S',
        'THR':'T', 'CYS':'C', 'TYR':'Y', 'ASN':'N', 'GLN':'Q',
        'ASP':'D', 'GLU':'E', 'LYS':'K', 'ARG':'R', 'HIS':'H'
                      }
    def __init__(self, aaSequence):
        self.aaSeq = aaSequence
        self.aaSeqLen = len(aaSequence)

    def cyclicSpectrum(self, aaSeq, cyclic=False):
        ''' TODO '''
        theoreticalSpectrum = [0]
        prefixMass = [0]

        # collect in a list the masses of first potential peptide (peptide prefix) substring based on length
        # for newAAseq in allPeptides:
        for aa in aaSeq:
            for testAA in self.aa2IntMass.keys():
                if testAA == aa:
                    # print(prefixMass)
                    # print(aa)
                    prefixMass.append(prefixMass[-1]+self.aa2IntMass[testAA])
        if cyclic:
            peptideMass = prefixMass[self.aaSeqLen] #collect the whole peptide mass

        # add to spectrum list the mass of everything else except the added peptide prefix
        for i in range(len(aaSeq)):
            for j in range(i+1,len(aaSeq)+1):
                # print(self.aaSeq[i:j+1])
                # print(i,j)
                theoreticalSpectrum.append(prefixMass[j]-prefixMass[i])
                # print(theoreticalSpectrum)
                if cyclic:
                    if i>0 and j<len(aaSeq):
                        testMass = peptideMass-(prefixMass[j]-prefixMass[i])
                        theoreticalSpectrum.append(testMass)

        # print(theoreticalSpectrum)
        theoreticalSpectrum.sort()
        return " ".join([str(x) for x in theoreticalSpectrum])
        # if len(theoreticalSpectrum) == (self.aaSeqLen*(self.aaSeqLen-1))+2:
        #     return " ".join([str(x) for x in theoreticalSpectrum])
        # else:
        #     return "ERROR OCCURED. SPECTRUM TOO LONG OR TOO SHORT FOR PEPTIDE SEQUENCE."

def main():
    for line in sys.stdin:
        myPeptide = CycloPeptide(line.rstrip())
        print(myPeptide.cyclicSpectrum(line.rstrip(),True))

if __name__ == "__main__": # if program is launched alone, this is true and is exececuted. if not, nothing is\
# executedf rom this program and instead objects and variables are made availableto the program that imports this.
    main();
    raise SystemExit
