/// Valuation.h
/// Shaun Harker
/// 2018-03-09
/// MIT LICENSE

#pragma once

#include "common.h"

std::function<Integer(Integer)>
construct_valuation ( std::shared_ptr<Complex> c, 
                      std::function<Integer(Integer)> top_cell_valuation ) {
  // Copy top_cell_valuation (with offset)
  std::vector<Integer> top_cell_valuation_;
  top_cell_valuation_.resize(c->size(c->dimension()));
  Integer num_nontop_cells_ = c->size() - c->size(c->dimension());
  for ( auto v : (*c)(c->dimension()) ) {
    top_cell_valuation_[v - num_nontop_cells_] = top_cell_valuation(v);
  }
  return [=](Integer x) { 
    Integer min_value = -1;
    for ( auto v : c->topstar(x) ) {
      auto new_val = top_cell_valuation_[v - num_nontop_cells_];
      if ( min_value == -1 ) {
        min_value = new_val;
      } else {
        min_value = std::min(min_value, new_val);
      }
    }
    return min_value;
  };
}

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline void
ValuationBinding(py::module &m) {
  m.def("construct_valuation", &construct_valuation);
}
