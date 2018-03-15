# InducedSubgraph.py
# Shaun Harker
# MIT LICENSE
# 2018-03-12

from pychomp.DirectedAcyclicGraph import *

def InducedSubgraph( G, predicate ):
    result = DirectedAcyclicGraph()
    S = set([v for v in G.vertices() if predicate(v)])
    for v in S:
        result.add_vertex(v)
    for v in S:
        for u in G.adjacencies(v):
            if u in S and u != v:
                result.add_edge(v,u)
    return result