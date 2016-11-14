#!/usr/bin/env python3
########################################################################
# File:treeBuilder.py
#  executable: treeBuilder.py <input >output
# Purpose: Build the Newick tree from tab delineated distance matrix.
#   stdout: Species tree in Newick format
#
# Author:       Akshar Lohith
# Group(sounding boards): Andrew Bailey
# History:      AL 11/10/2016 Created
#
#
#
########################################################################
import sys
class Node:
    '''Class to create node objects in a tree.
       Node objects are structures in a tree connected to other nodes, or not.

       Node class attributes:
        name: the name of the node
        outEdges: set of nodes that the node is a parent node.
        inEdges: set of nodes that have node as a child node.
        maxWeight: sum of the edge weights into the node.

       Node Class functions:
        printNode(): returns the node name.
        waysOut(): returns the set of nodes that call the node as parent
        addNodeConnection(direction, edge): update outEdges or inEdges
            sets with a edge object to connect.
        setNodeMaxWeight(): update the distance of the parent node.
        getMaxWeight(): return the distance of the parent node to the node.
    '''
    def __init__(self,name):
        self.name = name
        self.outEdges = list()
        self.inEdges = list()
        self.maxWeight = 0

    def printNode(self):
        '''Return the name of Node.'''
        return self.name

    def waysOut(self):
        '''Return the out edges of the Node.'''
        return self.outEdges

    def waysIn(self):
        '''Return the in edges of the Node.'''
        return self.inEdges

    def addEdge(self,direction,edge):
        '''Connect either an upstream or downstream node to the node.'''
        if direction == 'in':
            self.inEdges.append(edge)
        elif direction == 'out':
            self.outEdges.append(edge)

    def setNodeMaxWeight(self,weight):
        '''Update the sum of the edge weights into the node.'''
        self.maxWeight = weight

    def getMaxWeight(self):
        '''Return the sum of the edge weights into the node.'''
        return self.maxWeight

class Edge:
    '''Class to create Edge objects in a tree.
        Edge objects are structures in a tree connecting nodes to other nodes.

        Edge class attributes:
         outNode: the node object leading out of the edge.
         inNode: the node object leading into the edge.
         weight: the weight of the edge, given by problem statement.

        Edge Class functions:
         goesIntas(): returns the node object that is the input of the edge.
         goesOutas(): return the node object that is the output of the edge.
         giveWeight(): return the weight of the edge.
    '''
    def __init__(self,inNode,outNode, weight):
        self.inNode = inNode
        self.outNode = outNode
        self.weight = weight

        self.inNode.addEdge("out",self)
        self.outNode.addEdge('in',self)

    def goesIntas(self):
        '''Return the node that leads into the edge.'''
        return self.inNode

    def goesOutas(self):
        '''Return the node that leads out of the edge.'''
        return self.outNode

    def giveWeight(self):
        '''Return the weight of the edge.'''
        return self.weight

class Graph(object):
    """Class to create Graph objects. This class has functions that:
        createGraphObjects(edgeList): fill in the graph/DAG with node and edges.
        traceOutofSource(node): recursively fill the graph/DAG with weights, and
            sums of weights.
        tracebackFromSink(node): recursively back track the graph from sink to
            source node to obtain the highest scoring path.

        Class attributes:
        sourceNodeKey: key for the source node object
        sinkNodeKey: key for the sink node object
        nodeSet: dictionary containing node objects as values and node names as keys.
        listEdges: all possible edge objects
        graphMaxValue: the last node's maximum weight/value
        bestPathLength: the length/sum weight of the first highest scoring path.
        """
    def __init__(self, source, sink, edgeList):
        self.sourceNodeKey = source
        self.sinkNodeKey = sink
        self.nodeSet = dict()
        self.listEdges = list()
        self.graphMaxValue = 0
        self.bestPathLength = 0

        # fill in the graph/DAG
        self.createGraphObjects(edgeList)

        # the final path ends with the sink node, keep it in a list to add to in
        # traceback function.
        self.finalPath = [self.nodeSet[self.sinkNodeKey]]

    def createGraphObjects(self,edgeList):
        '''Given the adjacencyList and the source and sink node names, generate
            the graph/matrix/DAG containing node and edge objects with attributed
            weights.'''
        # start with the source and sink nodes
        self.nodeSet.update({self.sourceNodeKey:Node(self.sourceNodeKey),
            self.sinkNodeKey:Node(self.sinkNodeKey)})
        for edge in edgeList:
            nodeIN = edge[0]
            nodeOUT = edge[1]
            weight = edge[2]

            if nodeIN not in self.nodeSet.keys():
                self.nodeSet[nodeIN] = Node(nodeIN)
            if nodeOUT not in self.nodeSet.keys():
                self.nodeSet[nodeOUT] = Node(nodeOUT)

            self.listEdges.append(Edge(self.nodeSet[nodeIN],\
                    self.nodeSet[nodeOUT],int(weight)))

    def traceOutofSource(self,node):
        '''Starting from the Source node, transverse the graph/DAG and assign
            weights to the nodes based on the edges that lead up to that node.
            This is done recursively as every step possible in the graph is transversed.
            Credit Andrew Bailey for help with pseudocode that derived this code.'''
        if node == None:
            nodeList = [self.nodeSet[self.sourceNodeKey]]
        else:
            nodeList = node
        nextNodes = list()

        for node in nodeList:
            for edge in node.waysOut():
                nextNodes.append(edge.goesOutas())
                testWeight = node.getMaxWeight() + edge.giveWeight()

                if edge.goesOutas().getMaxWeight() < testWeight:
                    edge.goesOutas().setNodeMaxWeight(testWeight)
                    # regardless of step walked when the graphObj lands on the
                    # last node it will have the max score of the graph/matrix/DAG
                    self.graphMaxValue = testWeight

        if len(nextNodes) == 0:
            return self.graphMaxValue
        else:
            return self.traceOutofSource(nextNodes)

    def tracebackFromSink(self,node):
        '''Starting from the Sink node, traceback to the source of the graph/DAG
            the path that best led to the final node in the graph/DAG.
            This is done recursively as walking backward is inherently recursive
             (to me at least).
            Credit Andrew Bailey for help with pseudocode that derived this code.'''

        if node == None:
            testNode = self.nodeSet[self.sinkNodeKey]
        else:
            testNode = node

        for edge in testNode.waysIn():
            testWeight = edge.goesIntas().getMaxWeight()+edge.giveWeight()

            # immediateTracebackNode = edge.goesIntas().printNode()
            # print(testWeight, "is testWeight of",immediateTracebackNode+" and "+testNode.printNode(),"yo, TestNodeMax", testNode.getMaxWeight())

            if testNode.getMaxWeight()== testWeight:
                # the bestPathLength is the same as the penultimate node's weight,
                # plus the weight of the edge leading to the final node.
                # so set the best path length to the first node-edge pair that qualifies
                if node == None:
                    self.bestPathLength = testWeight

                # print("add {} to reverse path".format(edge.goesIntas().printNode()))
                self.finalPath.append(edge.goesIntas())
                if edge.goesIntas() != self.nodeSet[self.sourceNodeKey]:
                    self.tracebackFromSink(edge.goesIntas())
                else:
                    self.finalPath.reverse()
                    print(self.bestPathLength)
                    return '->'.join([node.printNode() for node in self.finalPath])
        return '->'.join([node.printNode() for node in self.finalPath])

def main():
    fileProvided = sys.stdin
    source = fileProvided.readline().rstrip()
    sink = fileProvided.readline().rstrip()
    edgeList = []
    for line in fileProvided:
        edgeList.append(line.rstrip().replace('->',':').split(':'))

    newGraph = Graph(source,sink,edgeList)
    largestValueInGraph = newGraph.traceOutofSource(None)
    longestPath = newGraph.tracebackFromSink(None)
    # print(longestPathLength)
    print(longestPath)

if __name__ == "__main__": # if program is launched alone, this is true and is exececuted. if not, nothing is\
# executedf rom this program and instead objects and variables are made availableto the program that imports this.
    main();
    # raise SystemExit
