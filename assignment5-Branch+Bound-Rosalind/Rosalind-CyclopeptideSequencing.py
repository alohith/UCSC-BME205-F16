#!/usr/bin/env python3
########################################################################
# File:Rosalind-CyclopeptideSequencing.py
#  executable: Rosalind-CyclopeptideSequencing.py <input >output
# Purpose: Find the aaMass chain of potential cyclic peptides for a given spectrum.
#   stdout: aaMass chain of potential cyclic peptides
#
# Author:       Akshar Lohith
# Group(sounding boards): Alison Tang, Andrew Bailey
# History:      AL 10/30/2016 Created
#
########################################################################
import sys
class CycloPeptide:
    """
       Class CycloPeptide handles peptide strings and spectrums for various outputs
       to solve Rosalind textbook track problems of the ba4 problem set stem.
       input:
           amino acid integer mass spectrum.

       class attributes:
           aa2IntMass: dictionary coding single amino acid string character to interger mass.
           aasingle2triple: dictionary converting aminoAcid notation from single-char to triple-char.
           aatriple2single: dictionary converting aminoAcid notation from triple-char to single-char.
           RNAcodonTable: dictionary converting 3-nucleotide RNA string to triple-char AA string.
           DNAcodonTable: dictionary converting 3-nucleotide DNA string to triple-char AA string.
           RNAcodon2AAsingle: dictionary converting 3-nucleotide RNA string to single-char AA string.
           DNAcodon2AAsingle: dictionary converting 3-nucleotide DNA string to single-char AA string.

       class functions:
            makeTheoreticalSpectrum(peptideSequence, boolean): boolean defaults to False, returning a
                spectrum for a linear peptide sequence.
            mass(peptide): returns the integer mass of a given peptide sequence.
            spectrumAApossible(): using the instantiated integer spectrum, returns
                the set of viable aminoAcids in the spectrum from aa2IntMass dictionary.
            cyclopeptideSequencing(): using a branch-bound algorithim, returns a set of
                peptide sequences that can be indicated by instantiated spectrum.
            aaMassChain(): takes the set of peptide sequences given by cyclopeptideSequencing()
                and converts them to aminoAcid integer mass chains. (ex. '113-128-186')
            peptideEncoding(nucleotideType,nucleotideSequence): returns to user the nucleotide
                sequences in all reading frames (including revComp) that codes for the aminoAcid
                given to object instantiation. (depreciated for cyclopeptideSequencing problem)

    """
    # Single aminoAcid to interger mass
    aa2IntMass = {'G':57,'A':71,'S':87,'P':97,'V':99,'T':101,'C':103,'I':113,'L':113,
              'N':114,'D':115,'K':128,'Q':128,'E':129,'M':131,'H':137,'F':147,
              'R':156,'Y':163,'W':186}

    aatriple2single = {
        'GLY':'G', 'ALA':'A', 'VAL':'V', 'LEU':'L', 'ILE':'I',
        'MET':'M', 'PHE':'F', 'TRP':'W', 'PRO':'P', 'SER':'S',
        'THR':'T', 'CYS':'C', 'TYR':'Y', 'ASN':'N', 'GLN':'Q',
        'ASP':'D', 'GLU':'E', 'LYS':'K', 'ARG':'R', 'HIS':'H',
        '---':'*'
                      }

    # As written, these are accessed as class attributes, for example:
    # CycloPeptide.aasingle2triple['A']
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

    RNAcodon2AAsingle = dict()
    DNAcodon2AAsingle = dict()

    for key, value in RNAcodonTable.items():
        RNAcodon2AAsingle.update({key:aatriple2single[value.upper()]})

    for key, value in DNAcodonTable.items():
        DNAcodon2AAsingle.update({key:aatriple2single[value.upper()]})

    def __init__(self, aaSequenceSpectrum):
        '''For this version of the CycloPeptide class an aminoAcid spectrum is
           passed to this constructor method. The spectrum string is converted
           to an interger list of containing each spectrum mass to be used by
           attribute functions.'''

        self.aaSpectrum = aaSequenceSpectrum.split(' ')
        self.intAASpectrum = [int(mass) for mass in self.aaSpectrum]

    def makeTheoreticalSpectrum(self, aaSeq, cyclic=False):
        '''Return a theoreticalSpectrum given an amino acid sequence. By default
           assumes that peptide sequence is linear. Setting cyclic to True returns
           theoretical spectrum for cyclic peptide. Coded based on Compeau and
           Pevzner pseudocode (2ndEd vol.1 pg 211-212)'''

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
            peptideMass = prefixMass[len(aaSeq)] #collect the whole peptide mass

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

    def mass(self,peptide):
        '''Return interger mass of peptide using aa2IntMass dictionary.'''

        mass = 0
        for aa in peptide:
            mass += self.aa2IntMass[aa]
        return mass

    def spectrumAApossible(self):
        """Return a set of possible aminoAcids given the interger spectrum and
           using the amino acid to interger mass dictionary."""

        possibleAA = set()
        for mass in self.intAASpectrum:
            for aa, intMass in self.aa2IntMass.items():
                if intMass == mass:
                    possibleAA.add(aa)
        return possibleAA

    def aaMassChain(self):
        """Return a set of peptides by interger mass chains from the
           cyclopeptideSequencing function."""

        massChainSet = set()
        for peptide in self.cyclopeptideSequencing():
            massChain = list()
            for aa in peptide:
                massChain.append(str(self.aa2IntMass[aa]))
            massChainSet.add('-'.join(massChain))
        return massChainSet

    def cyclopeptideSequencing(self):
        '''Return a set of cyclopeptides that can fit the spectrum. Code follows
           pseudocode given by Compeau and Pevzner. Assistance in completing this
           module given by Alison Tang, and Andrew Bailey'''

        # Initialize empty set containers and potential aa additions available given spectrum
        peptides = set()
        growingSeeds = self.spectrumAApossible()
        bestMatchingPeptides = set()
        peptideTreeNotBig = True

        # Seed the growable peptides with the potential aa additions available
        peptides.update(growingSeeds)
        # print(peptides)

        # Grow each branch while true
        while peptideTreeNotBig:
            # Have a temporary set containing the growable branches then grow the branches
            tempPeptides = peptides.copy()
            for peptide in tempPeptides:
                for aa in growingSeeds:
                    # print(aa)
                    peptides.add(peptide+aa)

            # time to trim branches that do not fit the spectrum
            # copy the recently grown peptides to discard if it doesn't fit the spectrum
            tempPeptides = peptides.copy()
            for peptide in peptides:
                # print(peptide)
                # print(self.mass(peptide))

                #generate a linear spectrum for each peptide
                testLinearSpectrum = set([int(mass) for mass in self.makeTheoreticalSpectrum(peptide).split(' ')])
                if self.mass(peptide) == max(self.intAASpectrum):
                    # only enter here if total mass of peptide is equal to max spectrum residue (the peptide mass)
                    testSpectrum = [int(mass) for mass in self.makeTheoreticalSpectrum(peptide, cyclic=True).split(' ')]
                    #check if the cyclospectrum of the peptide is same as the spectrum
                    # print(testSpectrum)
                    if testSpectrum==self.intAASpectrum:
                        # print(peptide)
                        # add this peptide branch to the viable cyclopeptide sequnces
                        bestMatchingPeptides.add(peptide)
                        # and get rid of the peptide from the temporary set so you don't run into it again
                        tempPeptides.discard(peptide)
                        # We are collecting viable cyclopeptide Sequences,
                        # you don't need to grow the branches anymore
                        peptideTreeNotBig =False

                #if the mass of the peptide is not equal to the max of the spectrum
                #check to see if the linear spectrum is even a viable branch in the spectrum
                elif not testLinearSpectrum.issubset(set(self.intAASpectrum)):
                    # trim the branch if it is not
                    tempPeptides.discard(peptide)

            # peptides growable are now the ones that have passed processing for ill fitting spectra
            peptides = tempPeptides

        return bestMatchingPeptides

def main():
    fileProvided = sys.stdin
    for line in fileProvided:
        spectrum = line.rstrip()
        myPeptide = CycloPeptide(spectrum)
        possiblePeptides = myPeptide.aaMassChain()
        for massList in possiblePeptides.copy():
            if len(possiblePeptides) !=1:
                print(massList, end=' ')
                possiblePeptides.discard(massList)
            else:
                print(massList)

if __name__ == "__main__": # if program is launched alone, this is true and is exececuted. if not, nothing is\
# executedf rom this program and instead objects and variables are made availableto the program that imports this.
    main();
    raise SystemExit
