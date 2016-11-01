#!/usr/bin/env python3
import sys
class KmerOperations:
    def __init__(self,kLen,kmerSet,genome):
        self.kLen = kLen
        if kmerSet == None:
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
            # print(nodes)
            # print(adjacencyList)
            #
            # if kmerSuffix not in adjacencyList.keys():
            #     adjacencyList[kmerSuffix] = list()
            #     # print(nodes[kmerPrefix].printNode(), nodes[kmerSuffix].printNode())
            #     # print(adjacencyList)
            #     adjacencyList[kmerPrefix].append(kmerSuffix)


            # self.kmerSet.discard(kmer)
            # # discard kmer as to not check against itself
            # for nextKmer in self.kmerSet:
            #     nextKmerPrefix = nextKmer[:-1]
            #     nextKmerSuffix = nextKmer[1:]
            #     if kmerSuffix == nextKmerPrefix:
            #         adjacencyList[kmer] = [nextKmer]
            # self.kmerSet.add(kmer)
        return adjacencyList

class Node:
    ''' docstring for Node.'''
    def __init__(self,name):
        self.name = name
        self.inEdges = list()
        self.outEdges = list()

    def printNode(self):
        '''Return the name of Node.'''
        return self.name

    def waysIn(self):
        '''Return the Edges that come into Node.'''
        return self.inEdges

    def waysOut(self):
        '''Return the Edges that leave from Node.'''
        return self.outEdges

    def addEdge(self,direction,edge):
        '''Connect an Edge in either the in or out directions.'''
        if direction == 'in':
            self.inEdges.append(edge)
        elif direction == 'out':
            self.outEdges.append(edge)
        else:
            raise Exception("Connecting Edge to Node with improper direction!\
                     Use direction = 'in' or 'out'.")

    def escapeRoutesLeft(self):
        '''Return True if out edges left in node, False if not.'''
        #Help from Andrew Bailey
        if False not in [edge.traversed for edge in self.outEdges]:
            return False
        else:
            return True

class Edge:
    ''' docstring for Edge. '''
    def __init__(self,inNode,outNode):
        #self.name = name # needed???
        self.inNode = inNode
        self.outNode = outNode
        self.traversed = False

        self.inNode.addEdge("out",self)
        self.outNode.addEdge("in",self)

    def goesIntas(self):
        '''Return out node for Edge.'''
        return self.inNode

    def goesOutas(self):
        '''Return in node for Edge.'''
        return self.outNode

    def beenTraversed(self):
        '''Change traversed boolean to True.'''
        self.traversed = True

    def traversedStatus(self):
        '''Return if Edge has been traversed.'''
        return self.traversed

class Graph:
    """docstring for Graph."""
    def __init__(self, streamIn):
        self.inFile = streamIn
        self.nodeList = dict()
        self.edges = list()
        self.makeNodesAndEdges(self.inFile)

    def makeNodesAndEdges(self,stream):
        '''  TODO '''
        for line in stream:
            startNode, endNode = (line.rstrip().split(' -> '))
            endNode = endNode.split(',')
            if startNode not in self.nodeList.keys():
                self.nodeList[startNode] = Node(startNode)
            for node in endNode:
                if node not in self.nodeList.keys():
                    self.nodeList[node] = Node(node)
                # print(nodeList[startNode].printNode(), nodeList[node].printNode())
                self.edges.append(Edge(self.nodeList[startNode],self.nodeList[node]))

    def walkThePath(self):
        ''' TODO '''
        import random
        startPathNode = random.choice([x for x in self.nodeList.values()])
        # print(startPathNode.printNode())
        runCycle = [startPathNode]
        while True:
            stepsAvailable = [edge for edge in startPathNode.waysOut() \
                                if not edge.traversedStatus()]
            step = random.choice(stepsAvailable)
            step.beenTraversed()
            # print(step.goesIntas().printNode(), step.goesOutas().printNode())
            runCycle.append(step.goesOutas())
            # print(len(runCycle))
            if step.goesOutas().escapeRoutesLeft():
                startPathNode = step.goesOutas()
            elif not step.goesOutas().escapeRoutesLeft():
                try:
                    newStartNode = random.choice([node for node in runCycle \
                                           if node.escapeRoutesLeft()])
                    # print(newStartNode.printNode())
                    runCycle = self.cleanPath(runCycle,newStartNode)
                    startPathNode = newStartNode

                except IndexError:
                    return runCycle
                    # print([node.printNode() for node in runCycle])
                    # return False

    def cleanPath(self, path, nextStart):
        ''' TODO '''
        clearPath = list()
        possiblePaths = list()

        # print("Time to Clean!!")
        newStartNodeIndex = len(path)-1-path[::-1].index(nextStart)
        possiblePaths.append(path[0:newStartNodeIndex+1])
        possiblePaths.append(path[newStartNodeIndex:-1])

        clearPath = possiblePaths[1]+possiblePaths[0]

        return clearPath

def main():

    aGraph = Graph(sys.stdin)
    theEulerianPath= aGraph.walkThePath()
    print('->'.join([node.printNode() for node in theEulerianPath]))

if __name__ == "__main__": # if program is launched alone, this is true and is exececuted. if not, nothing is\
# executedf rom this program and instead objects and variables are made availableto the program that imports this.
    main();
    raise SystemExit
