# TransitiveReduction.py
# Shaun Harker
# MIT LICENSE
# 2018-03-10

from pychomp.DirectedAcyclicGraph import *
from pychomp.TopologicalSort import *

def TransitiveReduction( G ):
    """ Return a new graph which is the transitive reduction of a DAG G """
    # Algorithm. Compute longest-path distance between each pair of vertices. 
    # Then, construct a graph consisting length-one longest paths.
    TS = TopologicalSort(G.vertices(), G.adjacencies)
    result = DirectedAcyclicGraph()
    for v in G.vertices():
        result.add_vertex(v)    
    for v in G.vertices():
        # Find single-source longest paths from v
        lp = { u : -1 for u in G.vertices() }
        lp[v] = 0
        for u in reversed(TS):
            val = lp[u]
            if val >= 0:
                for w in G.adjacencies(u):
                    if u != w:
                        lp[w] = max(val + 1, lp[w])
        for u in [ w for w in lp if lp[w] == 1 ]:
            result.add_edge(v, u)
    return result

