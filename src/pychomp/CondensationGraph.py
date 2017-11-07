### CondensationGraph.py
### MIT LICENSE 2016 Shaun Harker

from StronglyConnectedComponents import *
from DirectedAcyclicGraph import *
from collections import defaultdict

def CondensationGraph(vertices, adjacencies):
  components = StronglyConnectedComponents(vertices, adjacencies)
  mapping = defaultdict(int, { u : i for i, component in enumerate(components) for u in component })
  scc_dag = DirectedAcyclicGraph()
  for i, component in enumerate(components):
    #print("Examining SCC " + str(i))
    scc_dag.add_vertex(i)
    for u in component:
      #print("Examining adjacencies of vertex " + str(u) + " in SCC " + str(i) )
      for v in adjacencies(u):
        #print("Adjacency " + str(v) + " belongs to SCC " + str(mapping[v]) )
        if i != mapping[v]: scc_dag.add_edge(i,mapping[v])
  return scc_dag, mapping


