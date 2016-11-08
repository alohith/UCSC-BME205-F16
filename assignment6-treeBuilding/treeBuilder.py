#!/usr/bin/env python3
########################################################################
# File:treeBuilder.py
#  executable: treeBuilder.py <input >output
# Purpose: Build the Newick tree from tab delineated distance matrix.
#   stdout: Species tree in Newick format
#
# Author:       Akshar Lohith
# Group(sounding boards):
# History:      AL 11/04/2016 Created
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
        childLeafs: set of nodes that the node is a parent node.
            (terminal leaf nodes have no child nodes)
        parentBranch: set of nodes that have node as a child node.
            (root node has no parent nodes)
        parentDistance: distance to the parent node.

       Node Class functions:
        printNode(): returns the node name.
        childLeafsOut(): returns the set of nodes that call the node as parent
        addNodeConnection(direction, edge): update childLeafsOut or parentBranch
            sets with a node object to connect.
        updateParentDistance(): update the distance of the parent node.
        getParentDistance(): return the distance of the parent node to the node.
    '''
    def __init__(self,name):
        self.name = name
        self.childLeafs = set()
        self.parentBranch = set()
        self.parentDistance = 0.0

    def printNode(self):
        '''Return the name of Node.'''
        return self.name

    def childLeafsOut(self):
        '''Return the child leafs of the Node.'''
        return self.childLeafs

    def addNodeConnection(self,direction,edge):
        '''Connect either a child node or parent branch to the node.'''
        if direction == 'in':
            self.parentBranch.add(edge)
        elif direction == 'out':
            self.childLeafs.add(edge)

    def updateParentDistance(self,size):
        '''Update the distance of the node to the Parent Node.'''
        self.parentDistance = size

    def getParentDistance(self):
        '''Return the Parent node distance.'''
        return self.parentDistance

class TreeBuilder():
    """
        Class to create TreeBuilder objects, that build a tree of nodes.

        TreeBuilder class has attributes:
        distanceMatrix: the distance matrix detailing distance of nodes
            from each other. Initially provided upon creation.
        numSpecies: the number of species (nodes) in the tree.
        nodeSpecies: dictionary containing node objects with the node name as keys.
        nodesExhausted: dictionary to contain nodes and keys that have been merged
            into internal nodes. (just here to make sure that node objects
            are not destroyed after merging.)

        TreeBuilder class has functions:
        newick(node,firstNode): Node is the node to check/print, firstNode is
            boolean to tell if this is the root to start from.
        findTheRoot(): calls mergeNodes() until there is only 1 node left.
        makeQMatrix(): calculate a Q matrix from the size matrix.
        makeTreeSizeMatrix(): calculate the size matrix.
        calcBranchLength2Parent(posNode1, posNode2): uses the positions of the 2
            nodes that merge into the newNode to calculate their distance to the
            newNode.
        updateDistanceMatrix(posNode1,node1,posNode2,node2,newNode): Updates the
            distanceMatrix attribute to account for the merging of 2 nodes.
        findNodestoMerge(): Searches the Q matrix for nodes to merge.
        mergeNodes(): Merges 2 nodes from given by findNodestoMerge().
    """
    def __init__(self, distanceMatrix):
        self.distanceMatrix = distanceMatrix
        self.numSpecies = len(distanceMatrix)
        self.nodeSpecies = dict()
        self.nodesExhausted = dict()
        self.nodeSpecies.update({row[0]:Node(row[0]) for row in self.distanceMatrix})

    def newick(self,node,firstNode):
        '''Recursive function to return to format the tree in Newick format.
            format used: (A:0.1,B:0.2,(C:0.3,D:0.4):0.5);  distances and leaf names.
            Credit Ed Rice for pseudocode/temp code in TA office hours.'''
        if firstNode:
            # print(self.nodeSpecies.keys())
            node = self.nodeSpecies[node]

        childs = list(node.childLeafsOut())
        if len(childs) == 0:
            return node.printNode()
        else:
            return '('+self.newick(childs[0],False)+\
                ':{:.3f}'.format(childs[0].getParentDistance())+','\
                +self.newick(childs[1],False)+\
                ':{:.3f}'.format(childs[1].getParentDistance())+')'

    def findTheRoot(self):
        '''Functions to merge nodes from distanceMatrix until only root exists.
            Returns the key to root Node.'''
        rootNodeKey = None
        while self.numSpecies > 1:
            # print(self.numSpecies)
            rootNodeKey = self.mergeNodes()
        return rootNodeKey

    def makeTreeSizeMatrix(self):
        '''Return the size matrix of the nodes in the distanceMatrix by taking
            the sum of all columns in each row (associated with each node/species).'''
        sizeMatrix = list()
        for row in self.distanceMatrix:
            sizeMatrix.append([row[0],sum(row[1:])])
        return sizeMatrix

    def makeQMatrix(self):
        '''Return the score Q for each position in the distanceMatrix.
            Equation followed: Q[i,j]=(n-2)×d(i,j)-(S(i)+S(j))'''
        nodeSizes = self.makeTreeSizeMatrix()
        # print(nodeSizes, end = '\n\n\n')
        Qmatrix = list()

        for i, row in enumerate(self.distanceMatrix):
            # header = row[0]
            Qrow = list()
            for j, distance in enumerate(row[1:]):
                Qrow.append((self.numSpecies-2)*distance-(nodeSizes[i][1]+nodeSizes[j][1]))

            #keep the row Node name in the matrix for reference.
            Qrow.insert(0,row[0])
            Qmatrix.append(Qrow)

        return Qmatrix

    def findNodestoMerge(self):
        '''Using the QMatrix, return a list of mergeableNodes containing the
            first seen lowest QMatrix Node, the position of the matrix
            (the paired node), and the QMatrix score associated; and the
            final minimum QMatrix Score seen.'''
        Qmatrix2Eval = self.makeQMatrix()
        # print(Qmatrix2Eval, end= '\n\n')

        mergeableNodes = list()
        # set the initial minimum to be arbitrarily large
        rowMin = sys.maxsize
        for i, row in enumerate(Qmatrix2Eval):
            # go through each cell in the Q matrix
            # if the cell coordinates indicates that species is being compared
            # to itself, skip it
            # otherwise check to see if this Q score is lower than the last Q score
            # if it is, add to a list a tuple containg the species, the j-th cell
            # indicating the lowest Q score, and the Q score itself.
            for j, column in enumerate(row[1:]):
                # print(i,j,row,column)
                if i == j:
                    # print("passing")
                    pass
                elif column <= rowMin:
                    rowMin = column
                    mergeableNodes.append((row[0],j,column))
        [mergeableNodes.remove(potentialMerge) for potentialMerge in \
            mergeableNodes.copy() if potentialMerge[2] != rowMin]
        return mergeableNodes, rowMin

    def mergeNodes(self):
        '''Using the list of mergeableNodes and the lowest QMatrix value, merge
            the two nodes associated with the first appearing node with the
            lowest QMatrix value ONLY. Recalculate the distanceMatrix, remove
            the nodes merged to the nodesExhausted dictionary, and return the
            key of the node object to caller.'''
        mergeableNodes, lowest = self.findNodestoMerge()
        # print(mergeableNodes, lowest)
        newestNodeKey = ''
        for i in range(len(mergeableNodes)):
            checkNode1 = mergeableNodes[i]
            node1Min = checkNode1[2]
            for pos, row in enumerate(self.distanceMatrix):
                if row[0] == checkNode1[0]:
                    node1Pos = pos

            if node1Min == lowest:
                #make handles to refer to node names easier
                node1 = checkNode1[0]
                node2 = self.distanceMatrix[checkNode1[1]][0]
                node2Pos = checkNode1[1]
                newNodeName= node1 +'-'+ node2
                node1Obj = self.nodeSpecies[node1]
                node2Obj = self.nodeSpecies[node2]

                #create a new node that has a merged name, of the 2child nodes
                #and connect the child nodes to the new merged node with the distance
                newNodeObj = Node(newNodeName)
                newNodeObj.addNodeConnection('out',self.nodeSpecies[node1])
                newNodeObj.addNodeConnection('out',self.nodeSpecies[node2])
                node1ParentDist, node2ParentDist = self.calcBranchLength2Parent(\
                    node1Pos, node2Pos)

                # update the now childNodes with the parent node information
                self.nodeSpecies[node1].addNodeConnection('in',newNodeObj)
                self.nodeSpecies[node2].addNodeConnection('in',newNodeObj)
                self.nodeSpecies[node1].updateParentDistance(node1ParentDist)
                self.nodeSpecies[node2].updateParentDistance(node2ParentDist)

                #add the new merged node to the nodelist dictionary
                # and remove the nodes that resulted in the merged node
                # and update the distanceMatrix
                self.nodeSpecies.update({newNodeName:newNodeObj})
                for k in [node1,node2]:
                    self.nodesExhausted.update({k:self.nodeSpecies.pop(k)})
                self.updateDistanceMatrix(node2Pos,node2Obj,node1Pos,node1Obj,\
                        newNodeObj)
                newestNodeKey = newNodeName

                #stop after the first listed mergable node.
                break

        return newestNodeKey

    def calcBranchLength2Parent(self, pos1, pos2):
        '''Return the distance of each node being merged to new formed internal
            node.
            Follows equation: d(i,u)=(d(i,j)/2)+((1/(2(n-2)))*(S(i)-S(j)))'''
        sizeName1 = self.makeTreeSizeMatrix()[pos1][1]
        sizeName2 = self.makeTreeSizeMatrix()[pos2][1]

        # avoid DivisionZeroError by checking that the number of species nodes
        # to merge are greater than 2
        if self.numSpecies == 2:
            parentBranch2Name1 = (self.distanceMatrix[pos1][pos2+1]/2)+0
        else:
            parentBranch2Name1 = (self.distanceMatrix[pos1][pos2+1]/2)+\
                ((1/(2*(self.numSpecies-2)))*(sizeName1-sizeName2))

        return parentBranch2Name1, self.distanceMatrix[pos1][pos2+1]-parentBranch2Name1

    def updateDistanceMatrix(self, pos1,merged1, pos2,merged2, newNode):
        '''Update the distance Matrix for inclusion of the newNode, setting up
            the two nodes that make the newNode to be removed from future consideration.
            Follows equation: d(u,k)=(d(i,k)-d(i,u))/2+(d(j,k)-d(j,u))/2'''
        newDistanceMatrix = list()
        mergedRows = [newNode.printNode()]

        # print("nodeSpecies Keys:", self.nodeSpecies.keys())
        # print(merged1.printNode()+' Dist: '+str(merged1.getParentDistance()))
        # print(merged2.printNode()+' Dist: '+str(merged2.getParentDistance()))
        merged1ParentDist = merged1.getParentDistance()
        merged2ParentDist = merged2.getParentDistance()

        # first collect all the cells in the distance matrix that are not
        # affected by the 2 nodes merged in the original distanceMatrix format
        for i, row in enumerate(self.distanceMatrix):
            distanceRow = list()
            for j, column in enumerate(row[1:]):
                if (i == pos1 or i == pos2) or (j == pos1 or j == pos2):
                    pass
                else:
                    distanceRow.append(self.distanceMatrix[i][j+1])
            # add appropriately formated rows to new distanceMatrix
            if len(distanceRow) >0:
                newRow = [row[0]]
                [newRow.append(float(dist)) for dist in distanceRow]
                newDistanceMatrix.append(newRow)

        # do the actual new distance calculation for each node still in the matrix
        for j, column in enumerate(self.distanceMatrix[pos1][1:]):
            if j == pos1 or j == pos2:
                pass
            else:
                newDist =((self.distanceMatrix[pos1][j+1]-merged1ParentDist)/2)\
                    +((self.distanceMatrix[pos2][j+1]-merged2ParentDist)/2)
                mergedRows.append(newDist)

        # add the newly cacluated distances to each node to the end of row in matrix
        for j, newDist in enumerate(mergedRows[1:]):
            newDistanceMatrix[j].append(newDist)
            # print(newDistanceMatrix[j], end='\n\n')

        # add a 0.0 distance to end of the merged nodes distance from itself
        # and add to the end of the distance Matrix.
        # update the distance Matrix with the merged distance
        # reduce the species count in the distance matrix
        mergedRows.append(float(0))
        newDistanceMatrix.append(mergedRows)
        self.distanceMatrix = newDistanceMatrix
        self.numSpecies -= 1

def main():
    distanceMatrix = list()
    for line in sys.stdin:
        row = line.rstrip().split('\t')
        rowTemp = list()
        for column in row:
            try:
                rowTemp.append(int(column))
            except ValueError:
                rowTemp.append(column)
        distanceMatrix.append(rowTemp)
    # print(distanceMatrix, end='\n\n')
    newTree = TreeBuilder(distanceMatrix)
    treeRoot = newTree.findTheRoot()
    print(newTree.newick(treeRoot, True),end=';\n')

if __name__ == "__main__": # if program is launched alone, this is true and is exececuted. if not, nothing is\
# executedf rom this program and instead objects and variables are made availableto the program that imports this.
    main();
    # raise SystemExit
