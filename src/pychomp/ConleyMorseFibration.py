### ConleyMorseFibration.py
### MIT LICENSE 2016 Shaun Harker

from pychomp._chomp import *
from pychomp.CondensationGraph import *
from pychomp.StronglyConnectedComponents import *
from pychomp.DirectedAcyclicGraph import *
from pychomp.Poset import *

def ConleyMorseFibration(complex, discrete_flow):
  """
  Overview:
    Given a complex and a graph on its top dimensional cells,
    produce a Fibration such that the preimage of a down set
    is the collection of cells in the closure of all the 
    associated top cells
  Inputs:
    complex    : a complex
    flow_graph : an ordered pair (complex.cells(), edges)
                 where edge edge is an ordered pair of vertices
  Algorithm:
    Apply strongly connected components algorithm and determine
    reachability relation among the strong components to learn
    a poset. Associated to each poset vertex is a collection of
    top cells. 
  """

  # Step 1. Compute the poset of strongly connected components
  
  vertices = [ cell for cell in complex(complex.dimension())]
  (dag, mapping) = CondensationGraph(vertices, discrete_flow)
  #poset = Poset(dag)

  # Step 2. Extend the mapping from top-cells to all cells
  # Basic idea: since the component indexing furnishes a linear
  #             extension of the poset, we assign each cell to 
  #             the minimum indexed poset which contains a top cell
  #             it is incident.


  # for cell in reversed(range(0,len(complex))):
  #   current_value = mapping[cell]
  #   for bd_cell in complex.boundary({cell}):
  #     mapping[bd_cell] = min(mapping.get(bd_cell,current_value), current_value)

  def star(cell):
    result = set()
    stack = [ v for v in complex.coboundary({cell}) ]
    while stack:
      v = stack.pop()
      if v in result: continue
      result.add(v)
      for u in complex.coboundary({v}):
        stack.append(u)
    return result

  num_nontop_cells = complex.size() - complex.size(complex.dimension())

  valuation = lambda x : min([mapping[z] for z in star(x) if z >= num_nontop_cells])

  return dag, Fibration(complex, valuation) # lambda x : mapping[x])

  #return poset, chompy.Fibration(complex, lambda x : mapping[x])
