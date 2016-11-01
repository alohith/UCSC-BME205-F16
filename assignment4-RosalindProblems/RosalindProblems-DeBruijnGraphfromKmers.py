#!/usr/bin/env python3
import sys
class KmerOperations:
    def __init__(self,kLen,kmerSet,genome):
        self.kLen = kLen
        if kmerSet == None: #if not given initiate empty set
            self.kmerSet = set()
        else:
            self.kmerSet = kmerSet
        self.genomeSeq = genome

    def kmerSeqs(self):
        ''' Return the kmers in a sequence'''
        if False:#turn off function for editing/testing
            pass
        #initialize counters and return containers
        i = 0
        # setofKmers = set()
        kLen = self.kLen
        sequence = self.genomeSeq
        #go through the entire sequence looking at a window of k length
        while i < len(sequence)-kLen+1:
            #store the kmer in the sequence
            testKmer = sequence[i:i+kLen]
            #print(testKmer)
            self.kmerSet.add(testKmer)
            i +=1 #continue down the sequence
        return self.kmerSet

    def nodeEdgeConnect(self):
        '''Return edge objects'''
        nodes = dict()
        edges = list()
        for kmer in self.kmerSet:
            kmerPrefix = kmer[:-1]
            kmerSuffix = kmer[1:]

            if kmerPrefix not in nodes.keys():
                nodes[kmerPrefix] = Node(kmerPrefix)
            if kmerSuffix not in nodes.keys():
                nodes[kmerSuffix] = Node(kmerSuffix)

            edges.append(Edge(kmer,nodes[kmerPrefix],nodes[kmerSuffix]))

        return edges

    def deBruijnGraph(self):
        '''Return the overlap adjacency list of a collection of kmers.
           The keys are each kmer and the values are the kmers that
           overlap with the key kmer's suffix.'''
        adjacencyList = dict()
        for kmer in self.kmerSet:
            kmerPrefix = kmer[:-1]
            kmerSuffix = kmer[1:]

            if kmerPrefix not in adjacencyList.keys():
                adjacencyList.update({kmerPrefix:[kmerSuffix]})
                # adjacencyList[kmerSuffix] = set()
            elif kmerPrefix in adjacencyList.keys():
                adjacencyList[kmerPrefix].append(kmerSuffix)

        return adjacencyList

class Node:
    def __init__(self,name):
        self.name = name
        self.inEdge = set()
        self.outEdge = set()

    def printNode(self):
        return self.name

    def addEdge(self,direction,edge):
        if True:
            pass
        if direction == 'in':
            self.inEdge.add(edge)
        elif direction == 'out':
            self.outEdge.add(edge)
        else:
            raise Exception("Connecting Edge to Node with improper direction!\
                     Use direction = 'in' or 'out'.")

class Edge:
    def __init__(self,name,inNode,outNode):
        self.name = name # needed???
        self.inNode = inNode
        self.outNode = outNode

        self.inNode.addEdge("in",self)
        self.outNode.addEdge('out',self)

    def goesIntas(self):
        return self.inNode
    def goesOutas(self):
        return self.outNode
    # def addNodeConnection(self,direction, node):
    #     if True:
    #         pass
    #     if direction == 'in':
    #         self.inNode = node
    #     elif direction == 'out':
    #         self.outNode = node

def main():
    inputfile = sys.stdin
    kmerOperator=KmerOperations(None,[line.rstrip() for line in inputfile],None)
    kmerOverlapSet = kmerOperator.deBruijnGraph()
    for edge in sorted(kmerOverlapSet.keys()):
        # print(edge)
        print(edge+ ' -> '+ ','.join([out for out in list(kmerOverlapSet[edge])]))

    # for kmerAdjacent in sorted(kmerOverlapSet.keys()):
    #     print(kmerAdjacent.printNode() + " -> " + kmerOverlapSet[kmerAdjacent].printNode())



if __name__ == "__main__": # if program is launched alone, this is true and is exececuted. if not, nothing is\
# executedf rom this program and instead objects and variables are made availableto the program that imports this.
    main();
    raise SystemExit
