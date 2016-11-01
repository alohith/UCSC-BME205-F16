#!/usr/bin/env python3
def bestKmerSeq(sequence,kLen):
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

with open("rosalind_ba3a.txt",'r') as promptFile:
    firstLine = promptFile.readline().rstrip()
    sequenceLine = promptFile.readline().rstrip()
    with open("RosalindProblems-StringComposition-out.txt",'w') as writeFile:
        for kmer in sorted(bestKmerSeq(sequenceLine, int(firstLine))):
            writeFile.write(kmer+'\n')
