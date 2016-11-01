import math, random

class Scoring(object):
    """docstring for Scoring."""
    def __init__(self, pseudoCount, numSeq, kmerLen, kSeqList):
        self.pseudo = pseudoCount
        self.ksequences = kSeqList
        self.k = kmerLen
        self.numSeq = numSeq
        self.entropyScore = float(0)
        self.consensus = ''
        self.countMatrix = list()
        for i in range(0,self.k):
            self.countMatrix.append({'A':0,'G':0,'C':0,'T':0})

    def makeConsensus(self):
        '''docstring'''
        if False:
            pass
        profile = self.prProfile()
        for pos, column in enumerate(profile):
            baseInPos = ''
            scorebase = float(0)
            for base in profile[pos].keys():
                if profile[pos][base] > scorebase:
                    baseInPos = base
                    scorebase = profile[pos][base]
                # print((baseInPos,scorebase))
            self.consensus += baseInPos
            print(self.consensus)
        return self.consensus

    def prProfile(self):
        '''Return the probability scoring matrix that can be used to determine a consensus seq.'''
        if False:
            pass
        scoringMatrix = list()
        for column in self.makeCountMatrix():
            prColumn = {
                'A':float((self.pseudo+column['A'])/((4*self.pseudo)+self.numSeq)),
                'T':float((self.pseudo+column['T'])/((4*self.pseudo)+self.numSeq)),
                'G':float((self.pseudo+column['G'])/((4*self.pseudo)+self.numSeq)),
                'C':float((self.pseudo+column['C'])/((4*self.pseudo)+self.numSeq))}
            scoringMatrix.append(prColumn)
        return scoringMatrix

    def entropy(self):
        '''Return a float with the entropy score, as defined by the sum of the probabilities of each
        base in each position of the profile.'''
        if False:
            pass
        for column in self.prProfile():
            self.entropyScore += sum([x*math.log2(x) for x in column.values()])
        return self.entropyScore * -1

    def makeCountMatrix(self):
        '''Return a list of dictionaries, where each position in the list contains the count of
        each nucleotide.'''
        if False:
            pass
        for seq in self.ksequences:
            for i, base in enumerate(seq):
                self.countMatrix[i][base] += 1
        return self.countMatrix

def makeRandMotif():
    '''docstring'''
    if False:
        pass
    randomMotif = list()
    for seq in faSeqList:
        kSlice = random.randint(0,len(seq)-kLen)
        randomMotif.append(seq[kSlice:kSlice+kLen])
    return randomMotif

def bestKmerSeq(sequence, profile):
    ''' docstring'''
    if False:
        pass
    i = 0
    bestKmer = ''
    bestKmerScore = 0
    while i < len(sequence)-kLen+1:
        testKmer = sequence[i:i+kLen]
        score = 1
        for pos, base in enumerate(testKmer):
            score *= profile[pos][base]
        if score > bestKmerScore:
            bestKmer = testKmer
            bestKmerScore = score
        i +=1
    print((bestKmer,bestKmerScore))
    return (bestKmer, bestKmerScore)
def randomizedMotifSearch():
    ''' docstring'''
    if False:
        pass
    randomMotif = makeRandMotif()
    # make a motif (list of sequences of k length) randomly from the fasta sequences
    #create an object that has functions to make probability matrix, entropyScore, consensus sequence, and count matrix. using the 'bestMotif' as input
    evalRandomMotif = Scoring(pseudo,len(randomMotif),kLen,randomMotif)
    #generate the probability matrix (profile)
    prProfile1 = evalRandomMotif.prProfile()
    currentMotifEntropyScore = evalRandomMotif.entropy() #get bestMotif entropy score
    testMotifEntropyScore = 0
    bestMotif = randomMotif
    while True:
        testMotif = list()
        for seq in faSeqList:
            testMotif.append(bestKmerSeq(seq,prProfile1)[0])
        evalTestMotif = Scoring(pseudo,len(testMotif),kLen,testMotif)
        prProfile1 = evalTestMotif.prProfile()
        currentMotifEntropyScore = testMotifEntropyScore
        testMotifEntropyScore = evalTestMotif.entropy() #get newTestMotif entropy score
        if testMotifEntropyScore < currentMotifEntropyScore:
            bestMotif = testMotif
        else:
            evalBestMotif = Scoring(pseudo, len(bestMotif),kLen,bestMotif)
            return [bestMotif, evalBestMotif.entropy(), evalBestMotif.makeConsensus()]

def findBestMotif():
    '''docstring'''
    if False:
        pass
    # start a counter form 0 to the user inputted number of iterations
    iterCount = 0
    bestScore = 2*kLen
    bestMotifs = list()
    bestConsensus = ''
    #run user defined number of times
    while iterCount < usrIterations:
        testSet = randomizedMotifSearch()
        #test the two scores to see which is lower
        if testSet[1] < bestScore:
            bestScore = testSet[1]
            bestConsensus = testSet[2]
            bestMotif = testSet[0]
        #move to next iteration
        iterCount += 1
    return [bestMotif, bestScore, bestConsensus]


usrIterations = 10
kLen = 4
pseudo = 1
faSeqList = ['TTACCTTAAC',  'GATGTCTGTC', 'ACGGCGTTAG', 'CCCTAACGAG', 'CGTCAGAGGT']
