class Edge:
    def __init__(self,name,inA,outB):
        self.name = name
        self.inA = inA
        self.outB = outB
        self.fromNode.addOutEdge(self)
        self.toNode.addInEdge(self)

class Node:
    def __init__(self,name,):
        self.name = name
        self.incomingEdges = set()
        self.outgoingEdges = set()

# These objects may be overkill...
# want to build generic data structure that is the debruin graph


# 24 October 2016 TA/Tutor Stuff
nodes = dict()
edges = list()
for t in readEdges("edges.txt"):
    # update dictionary as you read your edges
    if t[0] not in nodes.keys():
        nodes[t[0]] = Node(t[0]) #value is the node object!! using the name
    if t[1] not in nodes.keys():
        nodes[t[1]] = Node(t[1])
    edges.append(Edge(nodes[t[0]],nodes[t[1]]) #edge knows about the nodes but not the other way
    nodes[t[0]].addOutEdge(edges[-1]) #this is nullified by last 2 lines in Edge init
    nodes[t[0]].addInEdge(edges[-1])
