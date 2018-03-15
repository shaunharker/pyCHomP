# TransitiveClosure.py
# Shaun Harker
# MIT LICENSE
# 2018-03-10

from pychomp.DirectedAcyclicGraph import *
from pychomp.TopologicalSort import *

def TransitiveClosure( G ):
    """ Return a new graph which is the transitive reduction of a DAG G """
    # Algorithm. Compute longest-path distance between each pair of vertices. 
    # Then, construct a graph consisting length-one longest paths.
    TS = TopologicalSort(G.vertices(), G.adjacencies)
    result = DirectedAcyclicGraph()
    for v in G.vertices():
        result.add_vertex(v)    
    for v in G.vertices():
        # Find vertices reachable from from v
        reachable = set()
        reachable.add(v)
        for u in reversed(TS):
            if u in reachable:
                for w in G.adjacencies(u):
                    reachable.add(w)
        for u in reachable:
            if u != v:
                result.add_edge(v, u)
    return result
