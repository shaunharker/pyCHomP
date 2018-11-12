/// MorseGradedComplex.h
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
#include "GradedComplex.h"

/// MorseGradedComplex
inline
std::shared_ptr<GradedComplex> 
MorseGradedComplex ( std::shared_ptr<GradedComplex> base_fibration, 
                 std::shared_ptr<MorseMatching> matching ) {

  std::shared_ptr<MorseComplex> complex ( new MorseComplex(base_fibration -> complex(), matching) );

  // Convert indices of cells to compute new fibration mapping (map from cell index to poset vertex number)
  std::vector<Integer> fibration_mapping(complex -> size());
  for ( auto x : *complex ) {
    Chain included = complex -> include ({x});
    fibration_mapping[x]= base_fibration -> value(*included.begin());
  }

  return std::shared_ptr<GradedComplex>( new GradedComplex(complex, [=](Integer x){return fibration_mapping[x];}));
}

/// MorseGradedComplex
inline
std::shared_ptr<GradedComplex> 
MorseGradedComplex ( std::shared_ptr<GradedComplex> base_fibration ) {
  std::shared_ptr<MorseMatching> matching ( MorseMatching::compute_matching(base_fibration) );
  return MorseGradedComplex (base_fibration, matching);
}

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline
void MorseGradedComplexBinding(py::module &m) {
  m.def("MorseGradedComplex", (std::shared_ptr<GradedComplex>(*)(std::shared_ptr<GradedComplex>,std::shared_ptr<MorseMatching>))&MorseGradedComplex);
  m.def("MorseGradedComplex", (std::shared_ptr<GradedComplex>(*)(std::shared_ptr<GradedComplex>))&MorseGradedComplex);
}
