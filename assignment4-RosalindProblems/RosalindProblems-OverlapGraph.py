#!/usr/bin/env python3
import sys
def kmerSeqs(sequence,kLen):
    ''' Return the kmers in a sequence'''
    if False:#turn off function for editing/testing
        pass
    #initialize counters and return containers
    i = 0
    setofKmers = set()
    #go through the entire sequence looking at a window of k length
    while i < len(sequence)-kLen+1:
        #store the kmer in the sequence
        testKmer = sequence[i:i+kLen]
        #print(testKmer)
        setofKmers.add(testKmer)
        i +=1 #continue down the sequence
    return setofKmers

def string2GenomePath(streamIn):
    ''' Return a reconstructed genome from its genome path.'''
    with streamIn as promptFile:
        finalGenome = promptFile.readline().rstrip()
        kmerLength = len(finalGenome)
        for line in promptFile:
            # print(line)
            # print(finalGenome)
            if line.rstrip()[:-1] == finalGenome[len(finalGenome)-kmerLength+1:]:
                finalGenome += line.rstrip()[-1]
            else:
                raise Exception("An issue arose from checking new kmer against end of genome!")
    return finalGenome

def overlap(kmerSet):
    '''Return the overlap adjacency list of a collection of kmers.
       The keys are each kmer and the values are the kmers that
       overlap with the key kmer's suffix.'''
    adjacencyList = dict()
    for kmer in kmerSet:
        kmerPrefix = kmer[:-1]
        kmerSuffix = kmer[1:]
        kmerSet.discard(kmer)
        for nextKmer in kmerSet:
            nextKmerPrefix = nextKmer[:-1]
            nextKmerSuffix = nextKmer[1:]
            if kmerSuffix == nextKmerPrefix:
                adjacencyList[kmer] = nextKmer
        kmerSet.add(kmer)
    return adjacencyList

kmerCollection = set()
for line in sys.stdin:
    kmerCollection.add(line.rstrip())
for adjacentPrefix, adjacentSuffix in overlap(kmerCollection).items():
    print(adjacentPrefix +' -> '+ adjacentSuffix)
