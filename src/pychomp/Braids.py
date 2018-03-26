### Braids.py
### MIT LICENSE 2016 Shaun Harker

from pychomp._chomp import *

from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

class BraidDiagram:
  def __init__(self, braid_skeleton):
    """
    Inputs:
      braid_skeleton : a list of lists such that 
         braid_skeleton[i][j] gives the value of strand i at position j
    Outputs:
      braid_specification: a tuple (n, m, x, pi) 
                           where n : number of strands 
                                 m : number of waypoints 
                                 x(i,j) means the value of strand i at position j
                                 pi is a permutation such that x(i,j+m) == x(pi(i),j)
    """
    self.n = len(braid_skeleton)
    self.m = len(braid_skeleton[0]) - 1
    self.permute_ = {}
    for i in range(0,self.n):
      for j in range(0,self.n):
        if braid_skeleton[i][self.m] == braid_skeleton[j][0]:
          self.permute_[i] = j
    self.braid_skeleton_ = braid_skeleton
    self.min_ = min([ min(braid_skeleton[i]) for i in range(0,len(braid_skeleton)) ])
    self.max_ = max([ max(braid_skeleton[i]) for i in range(0,len(braid_skeleton)) ])
    # convert
    self.thresholds = [ [float("-inf")] + sorted( [self.braid_skeleton_[i][j] for i in range(0,self.n)] ) + [float("inf")] for j in range(0,self.m) ]
    for j in range(0,self.m):
      self.thresholds[j][0] = self.thresholds[j][1] - 1.0
      self.thresholds[j][-1] = self.thresholds[j][-2] + 1.0


  def __call__(self, i, j):
    """
    Return the height of strand i at position j
    """
    return self.braid_skeleton_[i][j]

  def pi(self,i):
    """
    pi is a permutation such that x(i,j+m) == x(pi(i),j)
    """
    return self.permute_[i]

  def lap(self, coordinates):
    """
    Compute the lap number for a domain (given by a point in it)
    """
    #midpoints = [ sum(domain.bounds()[j]) / 2.0 for j in (list(range(0,self.m)) + [0]) ] 
    
    # on right fringe, give a big lap number
    if any( coord == self.n + 1 for coord in coordinates ):
      return 2 * self.n * self.m

    # otherwise, do the legitimate computation:
    midpoints = [ (self.thresholds[j][i+1] + self.thresholds[j][i])/2.0 for j,i in list(enumerate(coordinates)) + [(0,coordinates[0])]]
    return sum(self(i,j) <= midpoints[j] and self(i,j+1) >= midpoints[j+1] for j in range(0,self.m) for i in range(0,self.n))

  def draw(self, domain=None):
    x = np.arange(self.m+1)
    for i in range(0,self.n):
      plt.plot(x, [self(i,j) for j in range(0,self.m+1)])
    if domain:
      def f(x):
        if x[0] == -float("inf") and x[1] == -float("inf"):
          return self.min_ - 1.0
        if x[0] == float("inf") and x[1] == float("inf"):
          return self.max_ + 1.0
        if x[0] == -float("inf"):
          return x[1] - .5
        if x[1] == float("inf"):
          return x[0] + .5
        return (x[0] + x[1]) / 2.0
      strand = [ f(domain.bounds()[d]) for d in range(0, self.m) ]
      strand = strand + [strand[0]]
      plt.plot(x, strand, '--', color='b',)
    plt.show()

  def __repr__(self):
    self.draw()
    return "Braid Diagram"

def BraidComplex( braid_diagram ):
  """
  Overview:
    Given a specification for a "braids" dynamical system,
    return the associated cubical complex and flow graph.
  """

  # Unpack the input
  n = braid_diagram.n
  m = braid_diagram.m
  x = lambda i,j : braid_diagram(i,j)
  pi = lambda i : braid_diagram.pi(i)

  # Create the associated cubical complex

  thresholds = [ [float("-inf")] + sorted( [x(i,j) for i in range(0,n)] ) + [float("inf")] for j in range(0,m) ]
  # complex = CubicalComplex(CubicalGrid(thresholds))
  complex = CubicalComplex([ len(thresholds[j]) for j in range(0,m)])
  
  #lap = lambda x : braid_diagram.lap(complex.coordinates(x))

  lap_dict = {}
  def lap(x):
    if x not in lap_dict:
      lap_dict[x] = braid_diagram.lap(complex.coordinates(x))
    return lap_dict[x]

  # for x in complex(complex.dimension()):
  #   print( str(x) + " has coordinates " + str(complex.coordinates(x)) + " and lap number " + str(lap(x)))

  # Construct the domains
  # domains = [cell for cell in complex.cells() if cell.dimension() == m]
  # walls = [cell for cell in complex.cells() if cell.dimension() == m-1]

  domains = [cell for cell in complex(m)]
  walls = [cell for cell in complex(m-1)]

  # Construct the edge set
  edges = defaultdict(set)
  for wall in walls:
    # A wall can have 1 adjacent domain if it is off at infinity
    if len(complex.coboundary({wall})) == 1: continue
    # Otherwise, it has precisely 2 adjacent domains
    [u, v] = complex.coboundary({wall})
    if lap(u) <= lap(v):
      edges[v].add(u)
    if lap(u) >= lap(v):
      edges[u].add(v)

  # Identify collapsed strands
  collapsed_strands = [ i for i in range(0,n) if pi(i) == i ]

  #collapsed_vertices = [ CubicalCell([ [ x(i,j), x(i,j) ] for j in range(0,m) ]) for i in collapsed_strands ]
  collapsed_vertices_coords = [ [ braid_diagram.thresholds[j].index(braid_diagram(i,j)) for j in range(0,m) ] for i in collapsed_strands]
  collapsed_vertices = [ complex.cell_index(coordinates, 0) for coordinates in collapsed_vertices_coords ]

  # Connect all cubes in the star of any collapsed strand
  for v in collapsed_vertices:
    #print("collapsed vertex " + str(v) + " has coordinates " + str(complex.coordinates(v)) + " and shape " + str(complex.cell_shape(v)))
    #surrounding_walls = [ cell for cell in star(v) if cell.dimension() == m-1 ]
    surrounding_walls = [ cell for cell in complex.star({v}) if cell >= walls[0] and cell <= walls[-1] ]

    for wall in surrounding_walls:
      if len(complex.coboundary({wall})) == 1: continue
      [u, v] = complex.coboundary({wall})
      edges[u].add(v)
      edges[v].add(u)

  return (complex, lambda v : edges[v] )
