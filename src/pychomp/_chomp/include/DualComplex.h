/// DualComplex.h
/// Shaun Harker
/// 2018-03-12
/// MIT LICENSE

#pragma once

#include "common.h"

#include "Complex.h"

/// DualComplex
class DualComplex : public Complex {
public:

  DualComplex( std::shared_ptr<Complex> c ) : c_(c) {
    dim_ = c_ -> dimension();
    begin_.resize(dim_+2);
    Integer cumulative = 0;
    for ( Integer d = 0; d <= dim_; ++ d ) {
      begin_[d] = Iterator(cumulative);
      cumulative += c_ -> size(dim_ - d);
    }
    begin_[dim_ + 1] = c_ -> size();
  };

  /// column
  ///   Apply "callback" method to every element in ith column of
  ///   boundary matrix
  virtual void
  column ( Integer i, std::function<void(Integer)> const& callback) const final {
    auto transformed = [&](Integer x){ callback(size() - x - 1); };
    c_ -> row(size() - i - 1, transformed );
  };

  /// row
  ///   Apply "callback" method to every element in ith row of
  ///   boundary matrix
  virtual void
  row ( Integer i, std::function<void(Integer)> const& callback) const final {
    auto transformed = [&](Integer x){ callback(size() - x - 1); };
    c_ -> column(size() - i - 1, transformed );
  };

protected:
  std::shared_ptr<Complex> c_;
};

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline void
DualComplexBinding(py::module &m) {
  py::class_<DualComplex, std::shared_ptr<DualComplex>, Complex>(m, "DualComplex")
    .def(py::init<std::shared_ptr<Complex>>());
}
