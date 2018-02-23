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
#include "Fibration.h"

class MorseMatching {
public:
  /// mate
  virtual Integer mate ( Integer x ) const = 0;

  /// priority
  virtual Integer priority ( Integer x ) const = 0;

  /// compute_matching
  static
  std::shared_ptr<MorseMatching>
  compute_matching ( std::shared_ptr<Complex> complex );

  /// compute_matching
  static
  std::shared_ptr<MorseMatching>
  compute_matching ( std::shared_ptr<Fibration> fibration );
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
