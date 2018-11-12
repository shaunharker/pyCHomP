/// ConnectionMatrix.h
/// Shaun Harker
/// 2017-07-20
/// MIT LICENSE

#pragma once

#include <memory>

#include "Integer.h"
#include "Chain.h"
#include "Complex.h"
#include "MorseComplex.h"
#include "MorseMatching.h"
#include "GradedComplex.h"
#include "MorseGradedComplex.h"

/// ConnectionMatrix
inline
std::shared_ptr<GradedComplex> 
ConnectionMatrix ( std::shared_ptr<GradedComplex> base ) {
  std::shared_ptr<GradedComplex> next = base;
  do {
    base = next;
    next = MorseGradedComplex(base);
  } while ( next -> complex() -> size() != base -> complex() -> size() );
  return base;
}

/// ConnectionMatrix
inline
std::vector<std::shared_ptr<GradedComplex>>
ConnectionMatrixTower ( std::shared_ptr<GradedComplex> base ) {
  std::vector<std::shared_ptr<GradedComplex>> tower;
  std::shared_ptr<GradedComplex> next = base;
  std::shared_ptr<GradedComplex> last;
  do {
    tower.push_back(next);
    last = tower.back();
    next = MorseGradedComplex(last);
  } while ( next -> complex() -> size() != last -> complex() -> size() );
  return tower;
}

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline
void ConnectionMatrixBinding(py::module &m) {
  m.def("ConnectionMatrix", &ConnectionMatrix);
  m.def("ConnectionMatrixTower", &ConnectionMatrixTower);

}
