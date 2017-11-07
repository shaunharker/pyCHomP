/// Chain.h
/// Shaun Harker
/// 2017-07-19
/// MIT LICENSE

#pragma once

#include <unordered_set>
#include <iostream>

#include "Integer.h"

typedef std::unordered_set<Integer> Chain; 

inline void 
operator += ( Chain& lhs, Integer rhs ) { 
  if ( lhs.count(rhs) ) lhs.erase(rhs); else lhs.insert(rhs);
}

inline void 
operator += ( Chain& lhs, Chain const& rhs ) {
  for ( auto x : rhs ) lhs += x;
}

inline Chain
operator + ( Chain const& lhs, Chain const& rhs ) {
  Chain result;
  result += lhs;
  result += rhs;
  return result;
}

inline std::ostream &
operator << ( std::ostream & outstream, Chain const& print_me ) {
  outstream << "[";
  for ( auto x : print_me ) outstream << x << ",";
  outstream << "]";
  return outstream;
}

/// Python Bindings
// Note: using default STL wrapper for Chain

// #include <pybind11/pybind11.h>
// #include <pybind11/stl.h>
// namespace py = pybind11;
//
//
// void ChainBinding(py::module &m) {
//     py::class_<Chain>(m, "Chain")
//         .def(py::init<>())
//         .def(py::init<std::vector<uint64_t> const&>())
//         .def("__str__", )
//         .def("__repr__", )
//         .def("__add__", )
//         .def("__iadd__", )
//         .def("__sub__", )
//         .def("__isub__", )
//         .def("__iter__", )
//         .def("__len__", ).
// }
