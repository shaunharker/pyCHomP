# TopologicalSort.py
# MIT LICENSE 2016
# Shaun Harker

def TopologicalSort(vertices, adjacencies):
  """
  Topological Sort of a directed acyclic graph
  Example:
    vertices = [1,2,3,4,5]
    edges = [(1,2),(2,3),(2,4),(4,5),(3,5),(1,5)]
    adjacencies = lambda v :  [ j for (i,j) in edges if i == v ]
    print(TopologicalSort(vertices,adjacencies))
  """
  result = []
  preordered = set()
  postordered = set()
  def unvisited_children(u):
    return [ w for w in adjacencies(u) if w not in preordered ]
  for root in vertices:
    if root in preordered: continue
    stack = [root]
    def visit(u): 
      if u in preordered: 
        if u not in postordered:
          result.append(u)
          postordered.add(u)
      else:
        preordered.add(u)
        stack.extend([u] + unvisited_children(u))
    while stack: visit(stack.pop())
  return result
