/// OrderComplex.h
/// Shaun Harker
/// MIT LICENSE
/// 2018-03-12

#pragma once

#include "common.h"
#include "SimplicialComplex.h"

std::shared_ptr<SimplicialComplex>
OrderComplex ( std::shared_ptr<Complex> c ) {
  std::vector<Simplex> simplices;
  std::stack<Simplex> work_stack;
  for ( auto i : *c ) {
    work_stack.push({i});
    while ( not work_stack.empty() ) {
      Simplex s = work_stack.top();
      work_stack.pop();
      auto v = s.back();
      auto bd = c -> boundary({v});
      for ( auto u : bd ) {
        Simplex t = s;
        t.push_back(u);
        work_stack.push(t);
      }
      if ( bd.size() == 0 ) {
        simplices.push_back(s);
      }
    }
  }
  std::shared_ptr<SimplicialComplex> oc ( new SimplicialComplex(simplices) );
  return oc;
}

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline void
OrderComplexBinding(py::module &m) {
  m.def("OrderComplex", &OrderComplex);
}
