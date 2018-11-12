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
MorseGradedComplex ( std::shared_ptr<GradedComplex> base_graded_complex, 
                     std::shared_ptr<MorseMatching> matching ) {

  std::shared_ptr<MorseComplex> complex ( new MorseComplex(base_graded_complex -> complex(), matching) );

  // Convert indices of cells to compute new graded_complex mapping (map from cell index to poset vertex number)
  std::vector<Integer> graded_complex_mapping(complex -> size());
  for ( auto x : *complex ) {
    Chain included = complex -> include ({x});
    graded_complex_mapping[x]= base_graded_complex -> value(*included.begin());
  }

  return std::shared_ptr<GradedComplex>( new GradedComplex(complex, [=](Integer x){return graded_complex_mapping[x];}));
}

/// MorseGradedComplex
inline
std::shared_ptr<GradedComplex> 
MorseGradedComplex ( std::shared_ptr<GradedComplex> base_graded_complex ) {
  std::shared_ptr<MorseMatching> matching ( MorseMatching::compute_matching(base_graded_complex) );
  return MorseGradedComplex (base_graded_complex, matching);
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
