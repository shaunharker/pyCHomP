/// Valuation.h
/// Shaun Harker
/// 2018-03-09
/// MIT LICENSE

#pragma once

#include "common.h"

std::function<Integer(Integer)>
construct_grading ( std::shared_ptr<Complex> c, 
                    std::function<Integer(Integer)> top_cell_grading ) {
  // Copy top_cell_grading (with offset)
  std::vector<Integer> top_cell_grading_;
  top_cell_grading_.resize(c->size(c->dimension()));
  Integer num_nontop_cells_ = c->size() - c->size(c->dimension());
  for ( auto v : (*c)(c->dimension()) ) {
    top_cell_grading_[v - num_nontop_cells_] = top_cell_grading(v);
  }
  return [=](Integer x) { 
    // std::cout << "grading function\n";
    // std::cout << "top_cell_grading_.size() == " << top_cell_grading_.size() << "\n";
    // std::cout << "c -> size() == " << c -> size() << "\n";
    Integer min_value = -1;
    // for ( auto v : *c ) {
    //   std::cout << "c -> topstar(v) .size() == " << c ->topstar(v).size() << "\n";
    // }
    // std::cout << "grading of " << x << "\n";

    for ( auto v : c->topstar(x) ) {
      // std::cout << " top cell " << v << "\n";
      // std::cout << " minvalue = " << min_value << "\n";
      auto new_val = top_cell_grading_[v - num_nontop_cells_];
      // std::cout << "   new_val = " << new_val << "\n";
      if ( min_value == -1 ) {
        // std::cout << "    A\n";
        min_value = new_val;
      } else {
        // std::cout << "    B\n";
        min_value = std::min(min_value, new_val);
      }
    }
    // std::cout << "returning. minvalue = " << min_value << "\n";
    return min_value;
  };
}

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

namespace py = pybind11;

inline void
GradingBinding(py::module &m) {
  m.def("construct_grading", &construct_grading);
}
