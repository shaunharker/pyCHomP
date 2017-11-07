### StronglyConnectedComponents.py
### MIT LICENSE 2016 Shaun Harker

def StronglyConnectedComponents(vertices_input, adjacencies_input):
  """
  Overview:
    Compute the strongly connected components.
  Inputs:
    vertices_input    : a collection of vertices
    adjacencies_input : a function which takes a vertex and gives the adjacency list
  """

  # We translate the vertices to contiguous integers
  # via a list "vertices" and a dictionary "numbering"
  vertices = list(vertices_input)
  N = len(vertices)
  numbering = { v : i+1 for i, v in enumerate(vertices) }

  # We wrap the adjacencies list to use our contiguous integer indexing
  def Adjacencies(v):
    vertex = vertices[v-1]
    return [ numbering[u] for u in adjacencies_input(vertex)]

  # Data structures
  DFS = []
  S = []
  LOWLINK = [0] # dummy sentinel
  result = []
  index = {}
  explored = set()  #TODO, use array
  committed = set() #TODO, use array
  n = [0] # language caveat note: using list of one integer to capture variable and assign to it in "preorder" 
  # Preorder step 
  def preorder(v):
    if v in explored:
      return
    DFS.append(-v)
    explored.add(v)
    link = n[0]
    index[v] = link
    n[0] += 1
    for u in Adjacencies(v):
      if u not in explored:
        DFS.append(u)
      elif u not in committed:
        link = min(link, index[u])
    LOWLINK.append(link)
    S.append(v)

  # Postorder step 
  def postorder(v):
    lowlink = LOWLINK.pop()
    if lowlink == index[v]:
      SCC = []
      while v not in committed: #out of order? I had a do-while loop
        u = S.pop()
        SCC.append(u)
        committed.add(u)
      result.append(SCC)
    link = LOWLINK.pop()
    link = min(link, lowlink)
    LOWLINK.append(link)

  # Main routine
  DFS = list(range(1,N+1))
  while DFS:
    u = DFS.pop()
    if u > 0:
      preorder(u)
    else:
      postorder(-u)
  #Return result
  # Convert results to lists of vertices as originally given
  return [ [ vertices[v-1] for v in component ] for component in result ]

