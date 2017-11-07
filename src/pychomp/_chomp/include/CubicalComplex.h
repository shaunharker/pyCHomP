/// CubicalComplex.h
/// Shaun Harker 2016-12-15-2159
/// MIT LICENSE

#pragma once

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <numeric>
#include <vector>

#include <boost/functional/hash.hpp>

#include "Integer.h"
#include "Iterator.h"
#include "Chain.h"
#include "Complex.h"

/// CubicalComplex
///   Implements a trivial cubical complex with Z_2 coefficients
///   Methods:
///     ImageComplex : Initialize the complex with width N and height M
///     boundary : given a cell, returns an array of boundary cells
///     coboundary : given a cell, returns an array of coboundary cells
///     cells : return the set of cells in complex
class CubicalComplex : public Complex {
public:
  /// CubicalComplex
  ///   Default constructor
  CubicalComplex ( void ) {}

  /// CubicalComplex
  ///   Initialize the complex that is boxes[i] boxes across 
  ///   in dimensions d = 0, 1, ..., boxes.size() - 1
  ///   Note: The cubical complex does not have cells on the 
  ///         far right, so to have a "full" cubical 
  ///         complex as a subcomplex, pad with an extra box.
  CubicalComplex ( std::vector<Integer> const& boxes ) {
    assign ( boxes, std::vector<bool>(boxes.size(), false) );
  }

  /// CubicalComplex
  ///   Initialize the complex that is boxes[i] boxes across 
  ///   in dimensions d = 0, 1, ..., boxes.size() - 1
  ///   Note: The cubical complex does not have cells on the 
  ///         far right, so to have a "full" cubical 
  ///         complex as a subcomplex, pad with an extra box.
  CubicalComplex ( std::vector<Integer> const& boxes, 
                   std::vector<bool> const& periodic ) {
    assign ( boxes, periodic );
  }

  /// assign
  ///   Initialize the complex that is boxes[i] boxes across 
  ///   in dimensions d = 0, 1, ..., boxes.size() - 1
  void
  assign ( std::vector<Integer> const& boxes, std::vector<bool> const& periodic ) {
    // Get dimension
    Integer D = boxes.size();

    // Compute PV := [1, boxes[0], boxes[0]*boxes[1], ..., boxes[0]*boxes[1]*...*boxes[D-1]]
    auto & PV = place_values_;
    PV.resize ( D+1 );
    PV[0] = 1;
    std::partial_sum (boxes.begin(), boxes.end(), PV.begin() + 1, std::multiplies<Integer>()); 

    Integer M = 1L << D; // number of types
    Integer L = PV[D]; // Number of cells per shape/type
    Integer N = L * M; // cells per type * number of types

    // Set data
    boxes_ = boxes;
    dim_ = D;
    num_types_ = M;
    type_size_ = L;
    periodic_ = periodic;

    // Generate shapes and then sort them by dimension (which is bit popcount) to create types.
    // Implement a bijection between shapes and types via the following arrays:
    auto & ST = shape_from_type_;
    auto & TS = type_from_shape_;

    ST.resize ( M );
    TS.resize ( M );
    std::iota(ST.begin(), ST.end(), 0 );
    auto compare = [&](Integer x, Integer y) { return popcount_(x) < popcount_(y); };
    std::stable_sort(ST.begin(), ST.end(), compare);
    for ( Integer type = 0; type < M; ++ type) TS[ST[type]] = type; 

    // Set up iterator bounds for every dimension
    begin_ . resize ( dimension() + 2, N );
    for ( Integer type = 0, idx = 0; type < M; ++ type, idx += L ) {
      Integer shape = ST[type];
      Integer dim = popcount_(shape);
      begin_[dim] = Iterator(std::min(*begin_[dim], idx));
    }
  }

  /// boundary
  virtual Chain
  boundary ( Chain const& chain ) const final {
    Chain result;
    auto callback = [&](Integer bd_cell){result += bd_cell;};
    for ( auto x : chain ) column(x, callback);
    return result;
  }


  /// coboundary
  virtual Chain
  coboundary ( Chain const& chain ) const final {
    Chain result;
    auto callback = [&](Integer bd_cell){result += bd_cell;};
    for ( auto x : chain ) row(x, callback);
    return result;
  }

  /// column
  virtual void
  column ( Integer cell, std::function<void(Integer)> const& callback ) const final {
    Integer shape = cell_shape(cell);
    Integer type = TS() [ shape ];
    lldiv_t coordinate = {(int64_t)cell, 0}; // (quotient, remainder), see std::div
    for ( Integer d = 0, bit = 1; d < dimension(); ++ d, bit <<= 1L ) {
      // Determine dth coordinate
      coordinate = std::div(coordinate.quot, boxes()[d] ); 
      // If cell has no extent in this dimension, no boundaries.
      if ( not (shape & bit) ) continue;
      Integer offset_cell = cell + type_size() * ( TS() [ shape ^ bit ] - type );
      // Otherwise, the cell does have extent in this dimension.
      // It is always the case that such a cell has a boundary to the left.
      callback( offset_cell );
      // Check if there is a boundary to the right:
      if ( coordinate.rem + 1 < boxes()[d]) { 
        callback(offset_cell + PV()[d]);
      } else if ( periodic_[d] ) { // for periodic
        callback(offset_cell + PV()[d] - PV()[d+1]);     
      }
    }
  };

  /// row
  virtual void
  row ( Integer cell, std::function<void(Integer)> const& callback ) const final {
    Integer shape = cell_shape(cell);
    Integer type = TS() [ shape ];
    lldiv_t coordinate = {(int64_t)cell, 0};
    for ( Integer d = 0, bit = 1; d < dimension(); ++ d, bit <<= 1L ) {
      // Determine dth coordinate
      coordinate = std::div(coordinate.quot, boxes()[d] ); 
      // If cell has extent in this dimension, no coboundaries.
      if ( shape & bit ) continue;
      Integer offset_cell = cell + type_size() * ( TS() [ shape ^ bit ] - type );
      // Otherwise, the cell does not have extent in this dimension.
      // It is always the case that such a cell has a coboundary to the right.
      callback( offset_cell );
      // Check if there is a coboundary to the left:
      if ( coordinate.rem > 0 ) { 
        callback( offset_cell - PV()[d]);
      } else if ( periodic_[d] ) { // for periodic
        callback( offset_cell - PV()[d] + PV()[d+1]);
      }
    }
  };

  /// Features

  /// boxes
  ///   Number of boxes across in each dimension
  std::vector<Integer> const&
  boxes ( void ) const {
    return boxes_;
  }

  /// coordinates
  ///   Given a cell index, 
  ///   returns ( x_0, x_1, ..., x_{dim-1} )
  std::vector<Integer>
  coordinates ( Integer cell ) const {
    std::vector<Integer> result ( dimension() );
    for ( Integer d = 0; d < dimension(); ++ d ) {
      result[d] = cell % boxes()[d];
      cell /= boxes()[d];
    }
    return result;
  }

  /// shape_begin
  Iterator
  shape_begin ( Integer shape ) const {
    return Iterator(TS()[shape] * type_size());
  }

  /// shape_end
  Iterator
  shape_end ( Integer shape ) const {
    return Iterator(TS()[shape] * ( type_size() + 1));
  }

  /// cell_type
  ///   Give shape code
  ///   Interpretation: if ( shape & ( 1 << i ) ) { then the cell has extent in dimension i }
  Integer
  cell_type ( Integer cell ) const {
    return cell / type_size();
  }

  /// cell_shape
  ///   Give shape code
  ///   Interpretation: if ( shape & ( 1 << i ) ) { then the cell has extent in dimension i }
  Integer
  cell_shape ( Integer cell ) const {
    Integer shape = ST() [ cell_type(cell) ]; 
    return shape;
  }

  /// cell_index
  Integer
  cell_index ( std::vector<Integer> const& coordinates, 
               Integer shape ) {
    Integer cell = 0;
    for ( Integer d = dimension() - 1; d >= 0; -- d ) {
      cell *= boxes()[d];
      cell += coordinates[d];
    }
    cell += TS() [ shape ] * type_size();
    return cell;
  }

  /// periodic
  ///   Return vector reporting which dimensions
  ///   have periodic boundary conditions
  std::vector<bool> const&
  periodic ( void ) const {
    return periodic_;
  }

  /// cell_dim
  ///   Return dimension of cell
  Integer
  cell_dim ( Integer cell ) const {
    return popcount_(cell_shape(cell));
  }

  /// operator ==
  bool
  operator == ( CubicalComplex const& rhs ) const {
    if ( boxes() != rhs.boxes() ) return false;
    if ( periodic() != rhs.periodic() ) return false;
    return true;
  }

  /// operator <
  bool
  operator < ( CubicalComplex const& rhs ) const {
    if ( boxes() != rhs.boxes() ) {
      return std::lexicographical_compare(boxes().begin(), boxes().end(),
                                          rhs.boxes().begin(), rhs.boxes().end());
    } else {
      return std::lexicographical_compare(periodic().begin(), periodic().end(),
                                          rhs.periodic().begin(), rhs.periodic().end());
    }
  }

  /// operator <<
  friend std::ostream & operator << ( std::ostream & stream, CubicalComplex const& stream_me ) {
    stream << "CubicalComplex([";
    for ( auto x : stream_me.boxes() ) stream << x << ",";
    stream <<  "],["; 
    for ( auto x : stream_me.periodic() ) stream << ( x ? "T" : "F") << ",";
    return stream << "])";
  }

  /// print_cell (for debugging )
  void
  print_cell (Integer cell_index) const {
    std::cout << cell_index << "\n";
    std::cout << "  coordinates(" << cell_index << ") = [";
    for ( auto x : coordinates(cell_index) ) std::cout << x << ", "; std::cout << "]\n";
    std::cout << "  shape(" << cell_index << ") = " << cell_shape(cell_index) << "\n";
  };

private:

  std::vector<Integer> const&
  TS ( void ) const {
    return type_from_shape_;
  }

  std::vector<Integer> const&
  ST ( void ) const {
    return shape_from_type_;
  }

  std::vector<Integer> const&
  PV ( void ) const {
    return place_values_;
  }

  Integer
  type_size ( void ) const {
    return type_size_;
  }

  Integer
  popcount_ ( Integer x ) const {
    // http://lemire.me/blog/2016/05/23/the-surprising-cleverness-of-modern-compilers/
    Integer pcnt = 0; 
    while(x != 0) { x &= x - 1; ++pcnt; } 
    return pcnt;
  }
private:
  std::vector<Integer> boxes_;
  std::vector<Integer> place_values_;
  std::vector<Integer> shape_from_type_;
  std::vector<Integer> type_from_shape_;
  std::vector<bool> periodic_;
  Integer num_types_;
  Integer type_size_;
};

/// std::hash<CubicalComplex>
namespace std {
  template<> struct hash<CubicalComplex> {
    typedef CubicalComplex argument_type;
    typedef std::size_t result_type;
    result_type operator()(argument_type const& complex) const {
      using boost::hash_value;
      using boost::hash_combine;
      std::size_t seed = 0;
      for ( auto x : complex.boxes() ) {
        hash_combine(seed,hash_value(x));
      }     
      for ( auto x : complex.periodic() ) {
        hash_combine(seed,hash_value( x ? 1L : 0L ));
      }       
      return seed;
    }
  };
}

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline void
CubicalComplexBinding(py::module &m) {
  py::class_<CubicalComplex, std::shared_ptr<CubicalComplex>, Complex>(m, "CubicalComplex")
    .def(py::init<>())
    .def(py::init<std::vector<Integer> const&>())
    .def("coordinates", &CubicalComplex::coordinates)
    .def("cell_type", &CubicalComplex::cell_type)
    .def("cell_shape", &CubicalComplex::cell_shape)
    .def("cell_dim", &CubicalComplex::cell_dim)
    .def("cell_index", &CubicalComplex::cell_index);
}
