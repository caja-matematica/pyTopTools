"""
    
Metric for persistence diagrams
Author: [kel] 1/6/2013

"""

import hopcroft_karp as hk

### Added 2/18/2013
# Edited JJB: 5/7/13
def bdist ( d1, d2 ):
    """
    d1, d2 : list of tuple (birth, death)

    Returns distance and matching, else returns -1
    """

    if type( d1 ) == list:
        diagram_1 = d1[:]
    # if d1 is an array we have to convert so that extend() method
    # works to add generators
    else:
        diagram_1 = [ tuple( x ) for x in d1 ]
    if type( d2 ) == list:
        diagram_2 = d2[:]
    else:
        diagram_2 = [ tuple( x ) for x in d2 ]

    maxSize = len(diagram_1) + len(diagram_2)
    if maxSize == 0:
        return 0
    
    edges = []
    edges = prepareEdges(edges,diagram_1,diagram_2)
    
    graph = {}
    for i in xrange(2*maxSize):
        graph[i] = []
    
    matches = 0
    current_weight = 0
    first_not_added_edge = 0
    
    while matches < maxSize:
        while edges[ first_not_added_edge ].weight == current_weight and first_not_added_edge < len(edges)-1:
            graph [ edges[first_not_added_edge].vertex_1].append(edges[first_not_added_edge] . vertex_2 )
            graph [ edges[first_not_added_edge].vertex_2].append(edges[first_not_added_edge] . vertex_1 )
            first_not_added_edge += 1
        
        matching = hk.bipartiteMatch ( graph )[0]
        
        matches = len ( matching ) / 2
        optimalMatch = {}
        
        if matches == maxSize:
            for i in xrange(maxSize):
                optimalMatch [ diagram_1[i] ] = diagram_2 [ matching[i] - maxSize ]
            return current_weight, optimalMatch
            
        if first_not_added_edge == len (edges):
            print 'Not enough edges to find matching ? '
            return -1
        current_weight = edges[first_not_added_edge].weight
    return -1

class Edge:
    def __init__(self,v1,v2,w):
        self.vertex_1 = v1
        self.vertex_2 = v2
        self.weight = w
    def print_edge(self):
        print 'v1: ' + str(self.vertex_1) + ' v2: ' + str(self.vertex_2) + ' w: ' + str(self.weight)
        
def infDistanceOfTwoGenerators ( (b1,d1), (b2,d2) ):
    birth_diff = abs(b1-b2)
    death_diff = abs(d1-d2)
    if birth_diff > death_diff:
        return birth_diff
    return death_diff
def infDistanceOfGeneratorFromDiagonal ( (b,d) ):
    return abs(d-b)
    
def addGenerators ( diagram_1,diagram_2 ):
    newGens1 = []
    newGens2 = []
    for (b,d) in diagram_2:
        newGens1.append ( (b,b) )
    for (b,d) in diagram_1:
        newGens2.append ( (b,b) )
    diagram_1 . extend ( newGens1 )
    diagram_2 . extend ( newGens2 )

def prepareEdges( edges,diagram_1,diagram_2 ):
    
    maxSize = len(diagram_1)+len(diagram_2)
    #Connect diagonal points to each other
    for i in xrange(len(diagram_1), maxSize):
        for j in xrange(maxSize+len(diagram_2), 2*maxSize):
            edges.append ( Edge(i,j,0) )
    #Edges between non-diagonal points
    for i in xrange(len(diagram_1)):
        for j in xrange(maxSize,maxSize+len(diagram_2)):
            edges.append ( Edge(i,j,infDistanceOfTwoGenerators(diagram_1[i],diagram_2[j-maxSize]) ))
    #Edges between non-diagonal points and diagonal points
    for i in xrange(len(diagram_1)):
        edges.append( Edge(i,maxSize+len(diagram_2)+i,infDistanceOfGeneratorFromDiagonal(diagram_1[i])))
    
    for j in xrange(maxSize,maxSize+len(diagram_2)):
        edges.append( Edge(len(diagram_1)+(j-maxSize), j, infDistanceOfGeneratorFromDiagonal(diagram_2[j-maxSize])))
        
    edges = sorted(edges,key=lambda edge: edge.weight)
    
    #add diagonal generators to diagram for matching
    addGenerators(diagram_1,diagram_2)
    
    return edges

##### DEPRECATED ON 2/18/2013 ######

#from collections import deque
#import array
#### ----  BOTTLENECK DISTANCE  ---- ###
#
##HOPCROFT-KARP ALGORITHM#
##http://en.wikipedia.org/wiki/Hopcroft-Karp_algorithm
#
#def bfs(g1, pair_g1,pair_g2, adj, dist):
#    vertex_queue = deque([])
#    #Use 0 for index of NIL
#    for v in xrange(g1):
#        if pair_g1[v] == -1:
#            dist[v+1] = 0
#            vertex_queue.append(v)
#        else:
#            dist[v+1] = -1
#    dist[0] = -1
#    while len(vertex_queue) > 0:
#        v = vertex_queue.popleft()
#        for u in adj[v]:
#            if dist[ pair_g2[u]+1] == -1:
#                dist[pair_g2[u]+1] == dist[v+1] + 1
#                vertex_queue.append(pair_g2[u])
#    return (dist[0] != -1)
#    
#def dfs ( v, pair_g1, pair_g2, adj, dist ):
#    if v == -1:
#        return True
#    for u in adj[v]:
#        if dist[pair_g2[u]+1] == (dist[v+1] + 1):
#            if dfs(pair_g2[u], pair_g1,pair_g2,adj,dist):
#                pair_g2[u] = v
#                pair_g1[v] = u
#                return True
#    dist[v+1] = -1
#    return False
#
#def hopcroft_Karp (matching, maxSize,adj,dist,pair_g1,pair_g2):
#    while bfs(maxSize, pair_g1,pair_g2, adj, dist):
#        for vertex in xrange(maxSize):
#            if pair_g1[v] == -1:
#                if dfs(v, pair_g1,pair_g2, adj, dist):
#                    matching+=1
#    return matching
#        
#def bdistance ( d1, d2 ):
#    
#    diagram_1,diagram_2 = d1[:],d2[:]
#    maxSize = len(diagram_1) + len(diagram_2)
#    if maxSize == 0:
#        return 0
#    
#    edges = []
#    edges = prepareEdges(edges,diagram_1,diagram_2)
#    
#    adj = []
#    for i in xrange(maxSize):
#        adj.append([])
#    dist = array.array('i',(0,)*(maxSize+1))
#    pair_g1 = array.array('i',(-1,)*(2*maxSize))
#    pair_g2 = array.array('i',(-1,)*(2*maxSize))
#    
#    matching = 0
#    current_weight = 0
#    first_not_added_edge = 0
#    
#    while matching < maxSize:
#        while edges[ first_not_added_edge ].weight == current_weight and first_not_added_edge < len(edges)-1:
#            adj[edges[first_not_added_edge].vertex_1].append(edges[first_not_added_edge] . vertex_2 )
#            first_not_added_edge += 1
#        
#        matching = hopcroft_Karp(matching,maxSize,adj,dist,pair_g1,pair_g2)
#        
#        if matching == maxSize:
#            return current_weight
#            
#        if first_not_added_edge == len (edges):
#            print 'Not enough edges to find matching ? '
#            return -1
#        current_weight = edges[first_not_added_edge].weight
#    return -1
#




            
