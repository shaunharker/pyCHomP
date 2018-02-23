/// MorseFibration.h
/// Shaun Harker
/// 2017-07-20
/// MIT LICENSE

#pragma once

#include <memory>
#include <vector>

#include "Integer.h"
#include "Chain.h"
#include "Complex.h"
#include "MorseComplex.h"
#include "MorseMatching.h"
#include "Fibration.h"

/// MorseFibration
inline
std::shared_ptr<Fibration> 
MorseFibration ( std::shared_ptr<Fibration> base_fibration, 
                 std::shared_ptr<MorseMatching> matching ) {

  std::shared_ptr<MorseComplex> complex ( new MorseComplex(base_fibration -> complex(), matching) );

  // Convert indices of cells to compute new fibration mapping (map from cell index to poset vertex number)
  std::vector<Integer> fibration_mapping(complex -> size());
  for ( auto x : *complex ) {
    Chain included = complex -> include ({x});
    fibration_mapping[x]= base_fibration -> value(*included.begin());
  }

  return std::shared_ptr<Fibration>( new Fibration(complex, [=](Integer x){return fibration_mapping[x];}));
}

/// MorseFibration
inline
std::shared_ptr<Fibration> 
MorseFibration ( std::shared_ptr<Fibration> base_fibration ) {
  std::shared_ptr<MorseMatching> matching ( MorseMatching::compute_matching(base_fibration) );
  return MorseFibration (base_fibration, matching);
}

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline
void MorseFibrationBinding(py::module &m) {
  m.def("MorseFibration", (std::shared_ptr<Fibration>(*)(std::shared_ptr<Fibration>,std::shared_ptr<MorseMatching>))&MorseFibration);
  m.def("MorseFibration", (std::shared_ptr<Fibration>(*)(std::shared_ptr<Fibration>))&MorseFibration);
}
