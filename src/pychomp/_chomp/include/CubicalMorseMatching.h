/// CubicalMorseMatching.h
/// Shaun Harker
/// 2018-02-16
/// MIT LICENSE

#pragma once

#include <memory>
#include <unordered_set>
#include <vector>

#include "Integer.h"
#include "Chain.h"
#include "Complex.h"
#include "Fibration.h"
#include "MorseMatching.h"

class CubicalMorseMatching : public MorseMatching {
public:
  /// CubicalMorseMatching
  CubicalMorseMatching ( std::shared_ptr<CubicalComplex> complex_ptr ) : complex_(complex_ptr) {
    type_size_ = complex_ -> type_size();
    fibration_.reset(new Fibration(complex_, [](Integer i){return 0;}));
  }

  /// CubicalMorseMatching
  CubicalMorseMatching ( std::shared_ptr<Fibration> fibration_ptr ) : fibration_(fibration_ptr) {
    complex_ = std::dynamic_pointer_cast<CubicalComplex>(fibration_->complex());
    if ( not complex_ ) {
      throw std::invalid_argument("CubicalMorseMatching must be constructed with a Cubical Complex");
    }
    type_size_ = complex_ -> type_size();
  }

  /// mate
  Integer
  mate ( Integer x ) const { 
    return mate_(x, complex_ -> dimension());
  }

  /// priority
  Integer
  priority ( Integer x ) const { 
    return x % type_size_;
  }

private:
  uint64_t type_size_;
  std::shared_ptr<Fibration> fibration_;
  std::shared_ptr<CubicalComplex> complex_;

  // def mate(cell, D):
  // for d in range(0, D):
  //   if cell has extent in dimension d:
  //     left = leftboundary(cell, d)
  //     if value(left) == value(cell):
  //       if left == mate(left, d):
  //         return left
  //   else:
  //     right = rightcoboundary(cell, d)
  //     if value(right) == value(cell):
  //       if right == mate(right, d):
  //         return right
  //   return cell 
  // Note: should the complicated formulas (which are also found in CubicalComplex.h not be repeated here?
  Integer mate_ ( Integer cell, Integer D ) const {
    Integer shape = complex_ -> cell_shape(cell);
    Integer type = complex_ -> TS() [ shape ];
    for ( Integer d = 0, bit = 1; d < D; ++ d, bit <<= 1L  ) {
      if ( shape & bit ) { // if x has extent in dimension d
        Integer left_bd = cell + complex_ -> type_size() * ( complex_ -> TS() [ shape ^ bit ] - type );
        if ( fibration_ -> value(left_bd) == fibration_ -> value(cell) && left_bd == mate_(left_bd, d) ) 
          return left_bd;
      } else { // if x does not have extent in dimension d
        Integer right_cbd = cell + complex_ -> type_size() * ( complex_ -> TS() [ shape ^ bit ] - type );
        if ( fibration_ -> value(right_cbd) == fibration_ -> value(cell) && right_cbd == mate_(right_cbd, d) ) 
          return right_cbd;    
      }
    }
    return cell;
  }
};

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline void
CubicalMorseMatchingBinding(py::module &m) {
  py::class_<CubicalMorseMatching, std::shared_ptr<CubicalMorseMatching>>(m, "CubicalMorseMatching")
    .def(py::init<std::shared_ptr<CubicalComplex>>())
    .def(py::init<std::shared_ptr<Fibration>>())    
    .def("mate", &CubicalMorseMatching::mate)
    .def("priority", &CubicalMorseMatching::priority);
}
