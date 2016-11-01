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
    # firstLine = promptFile.readline().rstrip()
    # sequenceLine = promptFile.readline().rstrip()
    # with open("RosalindProblems-StringComposition-out.txt",'w') as writeFile:
    #     for kmer in sorted(bestKmerSeq(sequenceLine, int(firstLine))):
    #         writeFile.write(kmer+'\n')
print(string2GenomePath(sys.stdin))
