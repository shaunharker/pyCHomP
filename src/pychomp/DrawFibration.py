# DrawFibration.py

from collections import Counter
import graphviz

class DrawFibration:
    def __dir__(self):
        return list(self.__dict__.keys()) + dir(self._a)
    def __getattr__(self, attr):
        return getattr(self._a,attr)
    def __init__(self, a, poset):
        self._a = a
        self.poset = poset
        # Compute preimage
        self.preimage_ = {}
        for v in a.complex():
            val = a.value(v)
            if val not in self.preimage_:
                self.preimage_[val] = set()
            self.preimage_[val].add(v)
    def preimage(self, val):
        if val in self.preimage_:
            return self.preimage_[val]
        else:
            return set()
    def graphviz (self):
        """ Return a graphviz string describing the graph and its labels """
        gv = 'digraph {\n'
        indices = { v : str(k) for k,v in enumerate(self.poset.vertices())}
        counts = self._a.count()
        #print(counts)
        def vertex_label(v):
            if v in counts:
                return str(tuple(counts[v]))
            else:
                return " "
        for v in self.poset.vertices(): 
            gv += indices[v] + '[label="' + vertex_label(v) + ('", style=filled, fillcolor=cyan];\n' if self.preimage(v) else '"];\n')
        for v in self.poset.vertices(): 
            for u in self.poset.children(v):
                gv += indices[v] + ' -> ' + indices[u] + ';\n'
        return gv + '}\n'

    def _repr_svg_(self):
        return graphviz.Source(self.graphviz())._repr_svg_()
