### Poset.py
### MIT LICENSE 2016 Shaun Harker

import graphviz

class Poset:
  """
  Represents a poset
  """
  def __init__(self, graph):
    """  
    Create a Poset P from a DAG G such that x <= y in P iff there is a path from x to y in G 
    """
    self.vertices_ = set(graph.vertices())
    self.descendants_ = graph.transitive_closure()
    self.ancestors_ = self.descendants_.transpose()
    self.children_ = graph.transitive_reduction()
    self.parents_ = self.children_.transpose()

  def __iter__(self):
    """
    Allows for the semantics
      [v for v in poset]
    """
    return iter(self.vertices())

  def vertices(self):
    """ 
    Return the set of elements in the poset 
    """
    return self.vertices_
  
  def parents(self, v):
    """ 
    Return the immediate predecessors of v in the poset 
    """
    return self.parents_.adjacencies(v)
  
  def children(self, v):
    """ 
    Return the immediate successors of v in the poset 
    """
    return self.children_.adjacencies(v)
  
  def ancestors(self, v):
    """ 
    Return the set { u : u < v } 
    """
    return self.ancestors_.adjacencies(v)
  
  def descendants(self, v):
    """ 
    Return the set { u : v < u } 
    """
    return self.descendants_.adjacencies(v)
  
  def less(self, u, v):
    """ 
    Return True if u < v, False otherwise 
    """
    return u in self.ancestors(v)
  
  def maximal(self, subset):
    """ 
    Return the set of elements in "subset" which are maximal 
    """
    return frozenset({ u for u in subset if not any ( self.less(u,v) for v in subset ) })
  
  def _repr_svg_(self):
    """
    Return svg representation for visual display
    """
    return graphviz.Source(self.children_.graphviz())._repr_svg_()

def LatticeOfDownsets(poset):
  """ Generate from poset the Hasse diagram of the poset of down-sets of "poset" ordered by inclusion """
  lattice = Graph()
  recursion_stack = [poset.maximal(poset.vertices())]
  while recursion_stack:
    clique = recursion_stack.pop()
    for v in clique:
      parent_clique = poset.maximal(clique.difference([v]).union(poset.parents(v)))
      if parent_clique not in lattice.vertices():
        recursion_stack.append (parent_clique)
      lattice.add_edge (parent_clique,clique, str(v))
  return lattice