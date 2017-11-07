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
#include "Fibration.h"
#include "MorseFibration.h"

/// ConnectionMatrix
inline
std::shared_ptr<Fibration> 
ConnectionMatrix ( std::shared_ptr<Fibration> base ) {
  std::shared_ptr<Fibration> next = base;
  do {
    base = next;
    next = MorseFibration(base);
  } while ( next -> complex() -> size() != base -> complex() -> size() );
  return base;
}

/// ConnectionMatrix
inline
std::vector<std::shared_ptr<Fibration>>
ConnectionMatrixTower ( std::shared_ptr<Fibration> base ) {
  std::vector<std::shared_ptr<Fibration>> tower;
  std::shared_ptr<Fibration> next = base;
  std::shared_ptr<Fibration> last;
  do {
    tower.push_back(next);
    last = tower.back();
    next = MorseFibration(last);
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
