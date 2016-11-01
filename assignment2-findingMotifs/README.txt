README for randomizedMotifSearch.py:
usage: python3 randomizedMotifSearch.py [options] <input >output

This program will take from stdin a fasta file and return the consensus sequence based on the randomizedMotifSearch algorithim.

There are options to shuffle the sequences and return the consensus sequence and the entropy score of the shuffled sequences and
return the respective header and kmer for each of the sequences in the motif that contributed to the best consensus sequence.

Output from a run with non-default parameters:
$ python3 randomizedMotifSearch.py -p .2 -k 15 -r -i 2000 <catTestFiles
The Consensus of the Random Motif Search on the shuffled sequences is AAAACAAAAACCAAA -- with a score of 20.871925389717354
The Consensus of the Random Motif Search is ACTTAAAAACCCACA -- with a score of 9.598456296687235

Output from a run with default parameters:
$ python3 randomizedMotifSearch.py  <catTestFiles
The Consensus of the Random Motif Search is AAAAACTTAAAAA -- with a score of 12.760053874701878

Example output for help:
$ python3 randomizedMotifSearch.py -h
usage: randomizedMotifSearch.py [options] <input >output

This program will run a randomized motif search on a set of fasta records to
find the motif with the lowest entropy score, as defined by Shannon's entropy
formula

optional arguments:
  -h, --help            show this help message and exit
  -i NUMITERATE, --iterations NUMITERATE
                        Number of times to run the motif search before
                        returning a consensus sequence. Default is 1000
                        iterations.
  -p DEFINEDPSEUDOCOUNTS, --pseudoCount DEFINEDPSEUDOCOUNTS
                        Define a pseudo count to use for evaluation of the
                        kmers in the sequence. Default is 1 pseudo count.
  -k USRKLENGTH, --kmerLength USRKLENGTH
                        Set what kmer size to use to find the motif. This is a
                        required parameter.
  -r, --randomShuffle   Run the random motif search algorithim on shuffled
                        sequences.
  -m, --motifInRecord   Option to print the kmer contributing to consensus
                        from each sequence.