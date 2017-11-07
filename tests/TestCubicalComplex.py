import chomp
X = chomp.CubicalComplex([5,5])
print(len(X))
for cell in X:
  print(X.coordinates(cell))
  print(X.boundary(cell))
  print(X.coboundary(cell))
