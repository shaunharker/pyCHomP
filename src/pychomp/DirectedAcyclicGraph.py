### DirectedAcyclicGraph.py
### MIT LICENSE 2016 Shaun Harker

import subprocess, copy, json, graphviz, sys
from collections import defaultdict

from pychomp.TopologicalSort import *

# TODO: don't silently fail if given a non-DAG

class DirectedAcyclicGraph:
    """
    Represents a directed acyclic graph
    """
    def __init__(self):
        """ Initialize an empty graph object """
        self.vertices_ = set()
        self.adjacency_lists_ = {}
        self.vertex_labels_ = {}
        self.edge_labels_ = {}
    def add_vertex(self, v, label = ''):
        """ Add the vertex v to the graph and associate a label if one is given """
        if v in self.vertices_: return
        self.vertices_.add(v)
        self.adjacency_lists_[v] = set ()
        self.vertex_labels_[v] = label
        self.edge_labels_[v] = {}
    def add_edge(self, u, v, label = ''):
        """ Add the edge u -> v to the graph and associate a label if one is given """
        #print("Adding DAG edge (" + str(u) + ", " + str(v) + ")")
        self.add_vertex(u)
        self.add_vertex(v)
        self.adjacency_lists_[u].add(v)
        self.edge_labels_[u][v] = label
    def remove_edge(self, u, v):
        """ Remove the edge u -> v from the graph """
        self.adjacency_lists_[u].discard(v)
        self.edge_labels_[u].pop(v, None)
    def vertex_label(self, v):
        """ Return the label on the vertex v """
        return self.vertex_labels_[v]
    def get_vertex_from_label(self, label):
        """ Return the vertex v with label 'label'. Error if non-unique. """
        vertices = [ v for v in self.vertices_ if self.vertex_label(v) == label ]
        N = len(vertices)
        if N == 1:
            return vertices[0]
        elif N==0:
            return None
        elif N>1:
            raise ValueError("Non-unique vertex labels.")
    def edge_label(self, u, v):
        """ Return the label on the edge u -> v """
        return self.edge_labels_[u][v]
    def vertices(self):
        """ Return the set of vertices in the graph """
        return self.vertices_
    def edges(self):
        """ Return a complete list of directed edges (u,v) in the graph """
        return [(u,v) for u in self.vertices() for v in self.adjacencies(u)]
    def adjacencies(self, v):
        """ Return the set of adjacencies of v, i.e. { u : v -> u } """
        return self.adjacency_lists_[v]
    def clone(self):
        """ Return a copy of this graph """
        return copy.deepcopy(self)
    def transpose(self):
        """ Return a new graph with edge direction reversed. """
        G = DirectedAcyclicGraph ()
        for v in self.vertices(): G.add_vertex(v,self.vertex_label(v))
        for (u,v) in self.edges(): G.add_edge(v,u,self.edge_label(u,v))
        return G
    def transitive_closure(self):
        TS = TopologicalSort(self.vertices(), self.adjacencies)
        result = DirectedAcyclicGraph()
        for v in self.vertices():
            result.add_vertex(v)    
        for v in self.vertices():
            # Find vertices reachable from from v
            reachable = set()
            reachable.add(v)
            for u in reversed(TS):
                if u in reachable:
                    for w in self.adjacencies(u):
                        reachable.add(w)
            for u in reachable:
                if u != v:
                    result.add_edge(v, u)
        return result
        # """ Return a new graph which is the transitive closure """
        # G = self.clone ()
        # #print("Transitive closure: n = " + str(len(self.vertices())) )
        # for w in self.vertices():
        #     for u in self.vertices():
        #         for v in self.vertices():
        #             if w in G.adjacencies(u) and v in G.adjacencies(w):
        #                 G . add_edge(u,v)
        # return G
    def transitive_reduction(self):
        """ Return a new graph which is the transitive reduction """
        TS = TopologicalSort(self.vertices(), self.adjacencies)
        result = DirectedAcyclicGraph()
        for v in self.vertices():
            result.add_vertex(v)    
        for v in self.vertices():
            # Find single-source longest paths from v
            lp = { u : -1 for u in self.vertices() }
            lp[v] = 0
            for u in reversed(TS):
                val = lp[u]
                if val >= 0:
                    for w in self.adjacencies(u):
                        if u != w:
                            lp[w] = max(val + 1, lp[w])
            for u in [ w for w in lp if lp[w] == 1 ]:
                result.add_edge(v, u)
        return result

        # TC = self.transitive_closure ()
        # G = self.clone ()
        # for (u,v) in TC.edges():
        #     for w in TC.adjacencies(v):
        #         G.remove_edge(u,w)
        # return G
    def graphviz(self):
        """ Return a graphviz string describing the graph and its labels """
        gv = 'digraph {\n'
        indices = { v : str(k) for k,v in enumerate(self.vertices())}
        #for v in self.vertices(): gv += indices[v] + ';\n' #+ '[label="' + self.vertex_label(v) + '"];\n'
        for v in self.vertices(): gv += indices[v] + '[label="' + self.vertex_label(v) + '"];\n'
        for (u,v) in self.edges(): gv += indices[u] + ' -> ' + indices[v] + ' [label="' + self.edge_label(u,v) + '"];\n'
        return gv + '}\n'
    def _repr_svg_(self):
        return graphviz.Source(self.graphviz())._repr_svg_()
