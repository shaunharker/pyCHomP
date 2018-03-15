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
    //std::cout << "valuation function\n";
    //std::cout << "top_cell_valuation_.size() == " << top_cell_valuation_.size() << "\n";
    //std::cout << "c -> size() == " << c -> size() << "\n";
    Integer min_value = -1;
    // for ( auto v : *c ) {
    //   std::cout << "c -> topstar(v) .size() == " << c ->topstar(v).size() << "\n";
    // }
    //std::cout << "valuation of " << x << "\n";

    for ( auto v : c->topstar(x) ) {
      //std::cout << " top cell " << v << "\n";
      //std::cout << " minvalue = " << min_value << "\n";
      auto new_val = top_cell_valuation_[v - num_nontop_cells_];
      //std::cout << "   new_val = " << new_val << "\n";
      if ( min_value == -1 ) {
        //std::cout << "    A\n";
        min_value = new_val;
      } else {
        //std::cout << "    B\n";
        min_value = std::min(min_value, new_val);
      }
    }
    //std::cout << "returning. minvalue = " << min_value << "\n";
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
