/// MorseMatching.h
/// Shaun Harker
/// 2017-07-19
/// MIT LICENSE

#pragma once

#include <memory>
#include <unordered_set>
#include <vector>

#include "Integer.h"
#include "Chain.h"
#include "Complex.h"
#include "GradedComplex.h"

class MorseMatching {
public:
  typedef std::vector<Integer> BeginType; // To store location of first cell of each dim
  typedef std::vector<std::pair<Integer,Integer>> ReindexType; // To convert indexing

  /// mate
  virtual Integer mate ( Integer x ) const = 0;

  /// priority
  virtual Integer priority ( Integer x ) const = 0;

  /// critical_cells
  virtual std::pair<BeginType const&,ReindexType const&>
  critical_cells ( void ) const = 0;

  /// compute_matching
  static
  std::shared_ptr<MorseMatching>
  compute_matching ( std::shared_ptr<Complex> complex );

  /// compute_matching
  static
  std::shared_ptr<MorseMatching>
  compute_matching ( std::shared_ptr<GradedComplex> graded_complex );
};

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline void
MorseMatchingBinding(py::module &m) {
  py::class_<MorseMatching, std::shared_ptr<MorseMatching>>(m, "MorseMatching")  
    .def("mate", &MorseMatching::mate)
    .def("priority", &MorseMatching::priority);
}
