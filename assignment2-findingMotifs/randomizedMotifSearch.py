#!/usr/bin/env python3
########################################################################
# File:randomizedMotifSearch.py
#  executable: randomizedMotifSearch.py [options] <input >output
# Purpose: Find the lowest entropy motif using a randomizedMotifSearch
#   stderr: errors and status
#   stdout: the lowest entropy consensus sequence
#
# Author:       Akshar Lohith
# Group(sounding boards): Henry Gong, Brandon Saint-John, Bryan Thronlow
#                      Andrew Bailey IV
# History:      AL 10/05/2016 Created
#
#Program flow method
#   1) Implement randomizedMotifSearch as defined by the textbook on pg 93
#       a- take a number of sequences and create a motif of random kmers
#       b- create a profile (probability matrix) of this motifs and score it (entropyScore)
#       c- using this profile created, find a new motif with the most probable kmers of
#           of each sequence and score it.
#       d- evaluate for the lower entropy score.
#       e- if the new motif made has a lower entropy, set this motif (and it's profile)
#           as the best to test against the next time.
#       f- use this new best profile to find a new motif.
#       g- repeat from c until the new motif made is not of lower entropy.
#       h- return the motif with the lowest entropy and its consensus sequence.
#   2) Run part 1 for N (user defined or 1000 default) times collecting the
#       a final motif with the lowest entropy score and the consensus sequence
#       of all iterations.
#   3) Report this motif, and it's information to the user.
#Optionals completed:
#   1) Shuffle the sequences provided in sys.stdin and report the output of randomizedMotifSearch
#   2) Return the headers and kmers associated with them to the user
########################################################################
import sys
###########################################################################
# Scoring class
# here is where we evaluate motifs for the randomizedMotifSearch algorithim
###########################################################################
class Scoring(object):
    """Class to provide functions that help with the evaluation of
        motifs of k length.
        object instantiation:
        Scoring(<pseudoCount>,<length of motif>, <kmerLength>,<motif>)

        object attributes:
        pseudo: the provided pseudoCount to add to each element in countMatrix/profile
        ksequences: the motif (list of strings of k length) that is to be evaluated
        numSeq: the number of sequences in the motif. aka the number of fasta records user provides
        countMatrix: an empty countMatrix that is used to collect the count of each base in each position in the makeCountMatrix function

        methods:
        makeConsensus(): returns the consensus sequence of the motif given the probability scores of each base by position.
        prProfile(): returns the profile (probability matrix(probability of each base in each position)) of the motif.
        entropy(): returns the entropy score of the profile. calls prProfile().
        makeCountMatrix(): returns the countMatrix after adding the counts of each base in each position of the motif.

        Author: Akshar Lohith
        Date: 10 October 2016 (completed)
        """

    def __init__(self, pseudoCount, numSeq, kmerLen, kSeqList):
        '''constructor: saves attributes pseudoCount, numSeq, kmerLen, kSeqList, and creates empty countMatrix'''
        self.pseudo = pseudoCount
        self.ksequences = kSeqList
        self.k = kmerLen
        self.numSeq = numSeq
        self.countMatrix = list()
        #fill in the instantiated countMatrix with 0 for each base in each position
        for i in range(0,self.k):
            self.countMatrix.append({'A':0,'G':0,'C':0,'T':0})

    def makeConsensus(self):
        '''Returns the consensus sequence (string) based on the profile (prProfile()).'''
        # initialize the return containers
        consensus = ''
        profile = self.prProfile()

        #iterate over each position (k length) of the profile holding
        # both the dictionary in the position and the position index
        for pos, column in enumerate(profile):
            #initialize temporary container for the base and highest probability score
            baseInPos = ''
            scorebase = float(0)
            #evaluate the probability of each base in the index of the profile (the column of the matrix)
            # and add the base with the highest probability to the consensus sequence
            for base in profile[pos].keys():
                if profile[pos][base] > scorebase:
                    baseInPos = base
                    scorebase = profile[pos][base]
            consensus += baseInPos
        return consensus

    def prProfile(self):
        '''Return the probability scoring matrix that can be used to determine a consensus seq and entropy score.'''
        #initialize return container
        scoringMatrix = list()

        #calculate the probability of each base in the countMatrix by columns
        for column in self.makeCountMatrix():
            prColumn = {'A':float((self.pseudo+column['A'])/((4*self.pseudo)+self.numSeq)),
                        'T':float((self.pseudo+column['T'])/((4*self.pseudo)+self.numSeq)),
                        'G':float((self.pseudo+column['G'])/((4*self.pseudo)+self.numSeq)),
                        'C':float((self.pseudo+column['C'])/((4*self.pseudo)+self.numSeq))}
            scoringMatrix.append(prColumn)

        return scoringMatrix

    def entropy(self):
        '''Return a float with the entropy score, as defined by the sum of the probabilities of each
        base in each position of the profile.'''
        #initialize return container
        entropyScore = float(0)
        from math import log2 # to use log2 function

        #calculate the entropy Score by taking the sum of each base in each
        #column after changing to log base 2 probability
        for column in self.prProfile():
            entropyScore += sum([x*log2(x) for x in column.values()])

        #make the score positive, cause negative scores don't make sense
        return entropyScore * -1

    def makeCountMatrix(self):
        '''Return a list of dictionaries, where each position in the list contains the count of
        each nucleotide.'''
        #take each kmer in the motif
        #and add the count of each base to the matrix by position
        for seq in self.ksequences:
            for i, base in enumerate(seq):
                self.countMatrix[i][base] += 1
        return self.countMatrix

###########################################################################
# MotifSearch class
# here is where we create and operate the randomizedMotifSearch algorithim
###########################################################################
class MotifSearch(object):
    """Class that takes user's input to search for motifs based on algorithims in the textbook.

    Objection instantiation:
    MotifSearch(<fastaRecords>, <kmerLen>, <numIterate>, <pseudoCount>)

    object attributes:
    sequences: the list of fasta sequences to be used to search for motifs
    kLen: the user defined size of motif to search for.
    usrIterations: the number of times to run the Motif Search algorithim.
    pseudo: the user defined pseudoCount to keep as place holder in countMatrix.

    methods:
    makeRandMotif(): return a random motif of kmers from the fasta sequences in a list of sequences.
    bestKmerSeq(): return the most probable kmer in a sequence given a profile (scoring/probability matrix) and its probability in a tuple.
    randomizedMotifSearch(): return the bestMotif, its entropy score, and consensus sequence in a list.

    Author: Akshar Lohith
    Date: 10 October 2016 (completed)
    """

    def __init__(self,fastaRecords, kmerLen, numIterate, pseudoCount):
        '''constructor: saves the attributes fastaRecords, kmerLen, numIterate, pseudoCount'''
        self.sequences = fastaRecords
        self.kLen = kmerLen
        self.usrIterations = numIterate
        self.pseudo = pseudoCount

    def makeRandMotif(self):
        '''Returns a random motif of kmers from the fasta sequences provided.'''
        import random
        randomMotif = list()
        #go through each sequence one at a time and get an random index position
        # to slice out a kmer from the sequence, using the length of the sequence
        # and the desired kmer length
        for seq in self.sequences:
            kSlice = random.randint(0,len(seq)-self.kLen)
            randomMotif.append(seq[kSlice:kSlice+self.kLen])
        return randomMotif

    def bestKmerSeq(self, sequence, profile):
        ''' Return the FIRST highest probable kmer in a sequence given a probability matrix and its probability'''
        #initialize counters and return containers
        i = 0
        bestKmer = ''
        bestKmerScore = 0
        #go through the entire sequence looking at a window of k length
        #store the kmer in the sequence
        #go through the kmer base by base and get the kmer's probability based on the probability matrix
        while i < len(sequence)-self.kLen+1:

            testKmer = sequence[i:i+self.kLen]
            score = 1
            for pos, base in enumerate(testKmer):
                score *= profile[pos][base]

            #check the score of the kmer in the current window is better than the previous window
            if score > bestKmerScore:
                #store the kmer in the window and it's probability score if it is.
                bestKmer = testKmer
                bestKmerScore = score
            i +=1 #continue down the sequence

        return (bestKmer, bestKmerScore)

    def randomizedMotifSearch(self):
        ''' Implementation of the randomizedMotifSearch algorithim in textbook (pg 93). Return the motif with the lowest entropy score, its consensus sequence, and entropy score based on a randomly generated motif.'''
        # make a seeding motif (list of sequences of k length) randomly from the fasta sequences
        randomMotif = self.makeRandMotif()
        evalRandomMotif = Scoring(self.pseudo,len(randomMotif),self.kLen,randomMotif)

        #generate the seeding probability matrix (profile) to test and entropy score
        prProfile1 = evalRandomMotif.prProfile()
        currentMotifEntropyScore = evalRandomMotif.entropy()

        #initialize return containers
        testMotifEntropyScore = 0
        bestMotif = randomMotif
        while True: #run forever until the best motif is found
            #initialize a temporary list to hold the motif to check against
            testMotif = list()
            #get a motif of the most probable kmers given the currently best scoring profile
            for seq in self.sequences:
                testMotif.append(self.bestKmerSeq(seq,prProfile1)[0])
            evalTestMotif = Scoring(self.pseudo,len(testMotif),self.kLen,testMotif)

            #reset the currently best profile for next iteration and hold current best score
            prProfile1 = evalTestMotif.prProfile()
            currentMotifEntropyScore = testMotifEntropyScore
            testMotifEntropyScore = evalTestMotif.entropy()

            #test the currenlty best entropy score against the testMotif entropy score
            if testMotifEntropyScore < currentMotifEntropyScore:
                #if the test motif's entropy score is better (lower),
                #make the test motif the new best motif (with the best entropy score)
                bestMotif = testMotif
            else: #if not return the motif
                evalBestMotif = Scoring(self.pseudo, len(bestMotif),self.kLen,bestMotif)
                return [bestMotif, evalBestMotif.entropy(), evalBestMotif.makeConsensus()]

    def findBestMotif(self):
        '''Return the best motif, its entropyScore and consensus sequence after N (user defined) iterations'''
        # start a counter form 0 to the user inputted number of iterations
        iterCount = 0
        #set the 'best score' to the maximum score it could be
        bestScore = float(2*self.kLen)
        #initialize return containers
        bestMotifs = list()
        bestConsensus = ''

        #run the randomizedMotifSearch algorithim user defined number of times
        while iterCount < self.usrIterations:
            testSet = self.randomizedMotifSearch()

            #test the best and test scores to see which is of lower entropy
            # save the lower entropy score, motif, and consensus sequence
            if testSet[1] < bestScore:
                bestScore = testSet[1]
                bestConsensus = testSet[2]
                bestMotifs = testSet[0]

            iterCount += 1
        return [bestMotifs, bestScore, bestConsensus]

########################################################################
# FastAreader class
# here is where we parse out the fasta sequences from sys.stdin
########################################################################
class FastAreader :
    '''
    COPIED DIRECTLY FROM D. Bernick's BME160 LAB4 ASSIGNMENTS
    initial line modified to handle sys.stdin input

    Class to provide reading of a file containing one or more FASTA
    formatted sequences:
    object instantiation:
    FastAreader(<file name>):

    object attributes:
    fname: the initial file name

    methods:
    readFasta() : returns header and sequence as strings.
    Author: David Bernick
    Date: April 19, 2013
    '''
    def __init__ (self, fname):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def readFasta (self):
        '''
        using filename given in init, returns each included FastA record
        as 2 strings - header and sequence.
        whitespace is removed, no adjustment is made to sequence contents.
        The initial '>' is removed from the header.
        '''
        header = ''
        sequence = ''

        with self.fname as fileH:
            # initialize return containers
            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            # header is saved, get the rest of the sequence
            # up until the next header is found
            # then yield the results and wait for the next call.
            # next call will resume at the yield point
            # which is where we have the next header
            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()
        # final header and sequence will be seen with an end of file
        # with clause will terminate, so we do the final yield of the data
        yield header,sequence

########################################################################
# CommandLine class
# here is where we collect and parse out the commandline arguments
########################################################################
class CommandLine(object) :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    myCommandLine.args is a dictionary which includes each of the available command line arguments as
    myCommandLine.args['option']

    methods:
    do_usage_and_die()
    prints usage and help and terminates with an error.
    '''

    def __init__(self, inOpts=None) :
        '''
        CommandLine constructor.
        Implements a parser to interpret the command line argv string using argparse.
        '''
        import random
        import argparse
        self.parser = argparse.ArgumentParser(description = 'This program will run a randomized motif search on a set of fasta records to find the motif with the lowest entropy score, as defined by Shannon\'s entropy formula', add_help = True, #default is True
                                             prefix_chars = '-',
                                             usage = '%(prog)s [options] <input >output'
                                                  )
        self.parser.add_argument('-i', '--iterations', type=int, action = 'store', default=1000, dest='numIterate', help='Number of times to run the motif search before returning a consensus sequence. Default is 1000 iterations.')
        self.parser.add_argument('-p', '--pseudoCount', type=float, action = 'store', default=1, dest='definedPseudoCounts', help='Define a pseudo count to use for evaluation of the kmers in the sequence. Default is 1 pseudo count.')
        self.parser.add_argument('-k', '--kmerLength', type=int, action='store', default=13, dest='usrKlength',help='Set what kmer size to use to find the motif. This is a required parameter.')
        self.parser.add_argument('-r', '--randomShuffle', action='store_true', dest='shuffleOpt', help='Run the random motif search algorithim on shuffled sequences.')
        self.parser.add_argument('-m','--motifInRecord', action='store_true', dest='returnMotif', help='Option to print the kmer contributing to consensus from each sequence.')
        if inOpts is None :
            self.args = vars(self.parser.parse_args()) # parse the CommandLine options
        else :
            self.args = vars(self.parser.parse_args(inOpts)) # parse the input options


    def __del__ (self) :
        '''
        CommandLine destructor.
        '''
        # do something if needed to clean up before leaving
        pass

    def do_usage_and_die (self, str) :
        '''
        If a critical error is encountered, where it is suspected that the program is not being called with consistent parameters or data, this
        method will write out an error string (str), then terminate execution of the program.
        '''
        import sys
        print(str, file=sys.stderr)
        self.parser.print_usage()
        return 2

    def checkCommandLineOptions(self, inputSequences):
        '''Check the userProvided options and execute the MotifSearch algorithim.'''
        parseCmd = self.args
        sequencesFromPipe = inputSequences

        if parseCmd['returnMotif']:
            #run the randomizedMotifSearch algorithim user defined times, and collect the 'gems' (the information)
            #print out the header and the kmer sequence associated with that header in the motif together
            theGems = makeAndRunMotifSearch(sequencesFromPipe)
            print("The kmer sequences that made the consensus with lowest entropy are these sequences from each provided fasta record:")

            for header, motif in zip(headers,theGems[0]):
                print('>'+header)
                print(motif)

            #if the user specifies the shuffled option as well as the return the motif option do this..
            if parseCmd['shuffleOpt']:
                #make new container for shuffled sequences
                #shuffle function doesn't work on strings, so turn the sequence into a list
                #Shuffling happens in place

                seqShuffled = list()
                for seq in sequencesFromPipe:
                    tempSeq = list(seq)
                    random.shuffle(tempSeq)
                    seqShuffled.append(''.join(tempSeq))

                #call the motif search using the shuffled sequences as input
                #also call the motif search on the not-shuffled sequences for comparison
                theGems2 = makeAndRunMotifSearch(seqShuffled)
                #print a statement with the consensus and entropy score using a shuffled sequence input
                print("The Consensus of the Random Motif Search on the shuffled sequences is {} -- with a score of {}".format(theGems2[2],theGems2[1]))
            # print a statement with the consensus sequence and its entropy score
            print("The Consensus of the Random Motif Search is {} -- with a score of {}".\
                    format(theGems[2],theGems[1]))

        elif parseCmd['shuffleOpt']:

            seqShuffled = list()
            for seq in sequencesFromPipe:
                tempSeq = list(seq)
                random.shuffle(tempSeq)
                seqShuffled.append(''.join(tempSeq))

            #call the motif search using the shuffled sequences as input
            #also call the motif search on the not-shuffled sequences for comparison
            theGems2 = makeAndRunMotifSearch(seqShuffled)
            theGems = makeAndRunMotifSearch(sequencesFromPipe)
            #print a statement with the consensus and entropy score using a shuffled sequence input
            print("The Consensus of the Random Motif Search on the shuffled sequences is {} -- with a score of {}".\
                    format(theGems2[2],theGems2[1]))
            #print a statement with the consensus sequence and its entropy score
            print("The Consensus of the Random Motif Search is {} -- with a score of {}".\
                    format(theGems[2],theGems[1]))

        else:
            #if no options (print out motif, and shuffle the sequences) are raised, perform default fucntion of...
            #calling the motif search with the sequences as given and...
            theGems = makeAndRunMotifSearch(sequencesFromPipe)
            #printing a statement with the consensus sequence and its entropy score
            print("The Consensus of the Random Motif Search is {} -- with a score of {}".\
                    format(theGems[2],theGems[1]))

class Usage(Exception):
    '''
    Used to signal a Usage error, evoking a usage statement and eventual exit when raised.
    '''
    def __init__(self, msg):
        self.msg = msg

def makeAndRunMotifSearch(sequences,kmerLen=parseCmd['usrKlength'], \
        numIterate=parseCmd['numIterate'], pseudoCount=parseCmd['definedPseudoCounts']):
    '''Instantiate the MotifSearch object and return the output form the findBestMotif method'''

    findMyGems = MotifSearch(fastaRecords=sequences, kmerLen=kmerLen,\
                numIterate=numIterate, pseudoCount=pseudoCount)
                
    return findMyGems.findBestMotif()

########################################################################
# Main
# Here is the main program
########################################################################
def main(myCommandLine=None):
    '''
    Implements the Usage exception handler that can be raised from anywhere in process.
    '''
    if myCommandLine is None:
        myCommandLine = CommandLine()
    else :
        myCommandLine = CommandLine(['-h'])
    try:
        #initialize containers to parse the provided fasta records from stdin
        headers = list()
        sequencesFromPipe = list()
        # do the parsing
        for header, seq in FastAreader(sys.stdin).readFasta():
            headers.append(header)
            sequencesFromPipe.append(seq)

        myCommandLine.checkCommandLineOptions(sequencesFromPipe)

    except Usage as err:
       myCommandLine.do_usage_and_die(err.msg)

if __name__ == "__main__": # if program is launched alone, this is true and is exececuted. if not, nothing is\
# executedf rom this program and instead objects and variables are made availableto the program that imports this.
    main();
    raise SystemExit
