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
#Program flow method
########################################################################
import sys
class CycloPeptide:
    """docstring for CycloPeptide."""
    # Single aminoAcid to interger mass
    aa2IntMass = {'G':57,'A':71,'S':87,'P':97,'V':99,'T':101,'C':103,'I':113,'L':113,
              'N':114,'D':115,'K':128,'Q':128,'E':129,'M':131,'H':137,'F':147,
              'R':156,'Y':163,'W':186}
    aasingle2triple = {
        'G':'Gly', 'A':'Ala', 'V':'Val', 'L':'Leu', 'I':'Ile',
        'M':'Met', 'F':'Phe', 'W':'Trp', 'P':'Pro', 'S':'Ser',
        'T':'Thr', 'C':'Cys', 'Y':'Tyr', 'N':'Asn', 'Q':'Gln',
        'D':'Asp', 'E':'Glu', 'K':'Lys', 'R':'Arg', 'H':'His',
        '*':'---'
                      }
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
        # if len(theoreticalSpectrum) == (self.aaSeqLen*(self.aaSeqLen-1))+2:
        #     return " ".join([str(x) for x in theoreticalSpectrum])
        # else:
        #     return "ERROR OCCURED. SPECTRUM TOO LONG OR TOO SHORT FOR PEPTIDE SEQUENCE."

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

            # print(sorted(tempPeptides),'trimmedSet')
            # print(sorted(bestMatchingPeptides),'BestSet')
            # break

        return bestMatchingPeptides

    def peptideEncoding(self, nucleotideType,nucleotideSeq):
        '''Return the nucleotide sequnces in a DNA (or RNA) string that code for
           the given peptide sequnce.'''
        sequence = nucleotideSeq
        if nucleotideType == 'DNA':
            # initialize empty lists for pulling out aa sequences
            # and have the reverse complement sequence to process
            revCompSeq = rCompSeq(sequence)
            codonList = list()
            DNAStringConsumed = list()
            aaComplete = list()
            peptideSequenceFound = list()
            # go through each frame of sequence and pull out all possible aa sequences
            for frameCount in range(0,3):
                # go through the forward strand of the nucleotide sequence
                # print(frameCount, "+Frame {}".format(frameCount+1))
                for i in range(frameCount,len(sequence),3):
                    # print(sequence[i:i+3])
                    # print(i)
                    if sequence[i:i+3] in self.DNAcodon2AAsingle.keys():
                        aaComplete.append(self.DNAcodon2AAsingle[sequence[i:i+3]])
                        DNAStringConsumed.append(sequence[i:i+3])
                # print(DNAStringConsumed,aaComplete)
                codonList.append((DNAStringConsumed,aaComplete))

                # clear containers for processing the reverse strand
                DNAStringConsumed = list()
                aaComplete = list()
                # print(frameCount, "-Frame {}".format(frameCount+1))

            # do the same with the reverse complemented sequence at the same time
                for i in range(frameCount,len(revCompSeq),3):
                    # print(sequence[i:i+3])
                    # print(i)
                    if revCompSeq[i:i+3] in self.DNAcodon2AAsingle.keys():
                        aaComplete.append(self.DNAcodon2AAsingle[revCompSeq[i:i+3]])
                        DNAStringConsumed.append(revCompSeq[i:i+3])
                # print(DNAStringConsumed,aaComplete)
                codonList.append((DNAStringConsumed,aaComplete))
                DNAStringConsumed = list()
                aaComplete = list()
            # print(codonList)

            # go through the amino acid sequnce and find the aa sequecne given
            for j, element in enumerate(codonList):
                # print(element)
                for i in range(len(element[0])):
                    # print(len(element[0]))
                    # print(j, ''.join(element[1][i:i+self.aaSeqLen]))
                    if ''.join(element[1][i:i+self.aaSeqLen]) == self.aaSeq:
                        # join together the codons, revcomping the joined codons if on reverse strand frame
                        if j%2 ==1:
                            # print("".join(element[0][i:i+self.aaSeqLen]), rCompSeq("".join(element[0][i:i+self.aaSeqLen])),j,i)
                            peptideSequenceFound.append(rCompSeq("".join(element[0][i:i+self.aaSeqLen])))
                        else:
                            # print("".join(element[0][i:i+self.aaSeqLen]),j,i)
                            peptideSequenceFound.append("".join(element[0][i:i+self.aaSeqLen]))

            return peptideSequenceFound
        # do same as DNA but on RNA sequence:
        elif nucleotideType == 'RNA':
            # initialize empty lists for pulling out aa sequences
            revCompSeq = rCompSeq(sequence)
            codonList = list()
            RNAStringConsumed = list()
            aaComplete = list()
            peptideSequenceFound = list()
            # go through each frame of sequence and pull out all possible aa sequences
            for frameCount in range(0,3):
                # print(frameCount, "+Frame {}".format(frameCount+1))
                for i in range(frameCount,len(sequence),3):
                    # print(sequence[i:i+3])
                    # print(i)
                    if sequence[i:i+3] in self.RNAcodon2AAsingle.keys():
                        aaComplete.append(self.RNAcodon2AAsingle[sequence[i:i+3]])
                        DNAStringConsumed.append(sequence[i:i+3])
                # print(DNAStringConsumed,aaComplete)
                codonList.append((RNAStringConsumed,aaComplete))
                DNAStringConsumed = list()
                aaComplete = list()
                # print(frameCount, "-Frame {}".format(frameCount+1))

            # do the same with the reverse complemented sequence at the same time
                for i in range(frameCount,len(revCompSeq),3):
                    # print(sequence[i:i+3])
                    # print(i)
                    if revCompSeq[i:i+3] in self.RNAcodon2AAsingle.keys():
                        aaComplete.append(self.RNAcodon2AAsingle[revCompSeq[i:i+3]])
                        RNAStringConsumed.append(revCompSeq[i:i+3])
                # print(DNAStringConsumed,aaComplete)
                codonList.append((RNAStringConsumed,aaComplete))
                RNAStringConsumed = list()
                aaComplete = list()
            # print(codonList)

            for j, element in enumerate(codonList):
                # print(element)
                for i in range(len(element[0])):
                    # print(len(element[0]))
                    # print(j, ''.join(element[1][i:i+self.aaSeqLen]))
                    if ''.join(element[1][i:i+self.aaSeqLen]) == self.aaSeq:
                        if j%2 ==1:
                            # print("".join(element[0][i:i+self.aaSeqLen]), rCompSeq("".join(element[0][i:i+self.aaSeqLen])),j,i)
                            peptideSequenceFound.append(rCompSeq("".join(element[0][i:i+self.aaSeqLen])))
                        else:
                            # print("".join(element[0][i:i+self.aaSeqLen]),j,i)
                            peptideSequenceFound.append("".join(element[0][i:i+self.aaSeqLen]))
            return peptideSequenceFound

def rCompSeq(seq, nucType='DNA'):
    ''' Returns the reverse complement sequence of nucleotide string'''
    if nucType == 'DNA':
        bases = {'A':'T', 'T':'A', 'G':'C','C':'G'}
    elif nucType == 'RNA':
        bases = {'A':'U', 'U':'A', 'G':'C','C':'G'}
    return ''.join([x.replace(x,bases[x]) for x in list(seq)[::-1]])

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
