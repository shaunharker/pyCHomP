/// Fibration.h
/// Shaun Harker
/// 2017-07-20
/// MIT LICENSE

#pragma once

#include "Integer.h"
#include "Complex.h"
//#include "Poset.h"

class Fibration {
public:
  /// Fibration
  Fibration ( std::shared_ptr<Complex> c, 
              std::function<Integer(Integer)> v ) : complex_(c), value_(v) {}

  /// complex
  std::shared_ptr<Complex>
  complex ( void ) const {
    return complex_;
  }

  // /// poset
  // std::shared_ptr<Poset>
  // poset ( void ) const {
  //   return poset_;
  // }
  
  /// value
  Integer
  value ( Integer i) const {
    return value_(i);
  }

  /// count
  std::unordered_map<Integer,std::vector<Integer>>
  count ( void ) const {
    std::unordered_map<Integer,std::vector<Integer>> result;
    auto D = complex() -> dimension ();
    for ( Integer d = 0; d <= D; ++ d ) {
      for ( Integer idx : (*complex())(d) ) {
        auto v = value(idx);
        if ( result.count(v) == 0 ) result[v] = std::vector<Integer>(D+1);
        result[v][d] += 1;
      }
    }
    return result;
  }

private:
  std::shared_ptr<Complex> complex_;
  std::function<Integer(Integer)> value_;
};

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

namespace py = pybind11;

inline void
FibrationBinding(py::module &m) {
  py::class_<Fibration, std::shared_ptr<Fibration>>(m, "Fibration")
    .def(py::init<std::shared_ptr<Complex>,std::function<Integer(Integer)>>())
    .def("complex", &Fibration::complex)
    .def("value", &Fibration::value)
    .def("count", &Fibration::count);
}
