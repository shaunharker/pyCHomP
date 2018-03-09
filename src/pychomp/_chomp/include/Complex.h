/// Complex.h
/// Shaun Harker
/// 2017-07-18
/// MIT LICENSE

#pragma once

#include "common.h"

#include "Integer.h"
#include "Iterator.h"
#include "Chain.h"

/// Complex
class Complex {
public:

  /// virtual destructor
  virtual
  ~Complex ( void ) {}

  /// boundary
  virtual Chain
  boundary ( Chain const& chain ) const {
    Chain result;
    auto callback = [&](Integer bd_cell){result += bd_cell;};
    for ( auto x : chain ) column(x, callback);
    return result;
  }

  /// coboundary
  virtual Chain
  coboundary ( Chain const& chain ) const {
    Chain result;
    auto callback = [&](Integer bd_cell){result += bd_cell;};
    for ( auto x : chain ) row(x, callback);
    return result;  
  }

  /// star
  virtual std::unordered_set<Integer>
  star ( Integer cell ) const {
    std::unordered_set<Integer> result;
    std::stack<Integer> work_stack;
    work_stack.push(cell);
    while ( not work_stack.empty() ) {
      auto v = work_stack.top();
      work_stack.pop();
      if ( result.count(v) ) continue;
      result.insert(v);
      for ( auto u : coboundary({v}) ) {
        work_stack.push(u);
      }
    }
    return result;
  }

  /// topstar
  ///   return top dimensional cells in star
  virtual std::vector<Integer>
  topstar ( Integer cell ) const {
    Integer N = size() - size(dimension());
    std::vector<Integer> result;
    for ( auto v : star(cell) ) {
      if ( v >= N ) result.push_back(v);
    }
    return result;
  }

  /// column
  ///   Apply "callback" method to every element in ith column of
  ///   boundary matrix
  virtual void
  column ( Integer i, std::function<void(Integer)> const& callback) const {};

  /// row
  ///   Apply "callback" method to every element in ith row of
  ///   boundary matrix
  virtual void
  row ( Integer i, std::function<void(Integer)> const& callback) const {};
  
  /// dimension
  Integer 
  dimension ( void ) const {
    return dim_;
  }

  /// begin
  Iterator
  begin ( void ) const {
    return begin_[0];
  }

  /// end
  Iterator
  end ( void ) const  {
    return begin_.back();
  }

  /// size
  Integer
  size ( void ) const {
    return * begin_.back();
  }

  /// size
  Integer
  size ( Integer d ) const {
    if ( d < 0  || d > dim_ ) return 0;
    return begin_[d+1] - begin_[d];
  }

  /// operator ()
  Range
  operator () ( Integer dim ) const {
    return Range(begin_[dim], begin_[dim+1]);
  }

  /// count
  std::vector<Integer>
  count ( void ) const {
    std::vector<Integer> result;
    auto D = dimension ();
    for ( Integer d = 0; d <= D; ++ d ) result.push_back(size(d));
    return result;
  }

protected:
  Integer dim_;
  std::vector<Iterator> begin_; // begin_by_dim_[D+1] == size_;
};

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline void
ComplexBinding(py::module &m) {
  py::class_<Complex, std::shared_ptr<Complex>>(m, "Complex")
    .def("dimension", &Complex::dimension)
    .def("boundary", &Complex::boundary)
    .def("coboundary", &Complex::coboundary)
    .def("column", &Complex::column)
    .def("row", &Complex::row)
    .def("star", &Complex::star)
    .def("topstar", &Complex::topstar)
    .def("__iter__", [](Complex const& v) {
       return py::make_iterator(v.begin(), v.end());
    }, py::keep_alive<0, 1>())
    .def("__call__", [](Complex const& v, Integer d) {
       return py::make_iterator(v(d).begin(), v(d).end());
    }, py::keep_alive<0, 1>())
    .def("__len__", (Integer(Complex::*)(void)const)&Complex::size)
    .def("size", (Integer(Complex::*)(void)const)&Complex::size)
    .def("size", (Integer(Complex::*)(Integer)const)&Complex::size)
    .def("count", &Complex::count);
}
