/// CubicalComplex.h
/// Shaun Harker 2016-12-15-2159
/// MIT LICENSE

#pragma once

#include "common.h"

#include "Integer.h"
#include "Iterator.h"
#include "Chain.h"
#include "Complex.h"

// TODO:
//  get rid of interfaces like boxes()[d] in favor of boxes(d)

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
    assign ( boxes );
  }

  /// assign
  ///   Initialize the complex that is boxes[i] boxes across 
  ///   in dimensions d = 0, 1, ..., boxes.size() - 1
  void
  assign ( std::vector<Integer> const& boxes ) {
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

    // Generate shapes and then sort them by dimension (which is bit popcount) to create types.
    // for d = 4: (0D) 0, (1D) 1, 2, 4, 8, (2D) 3, 5, 6, 9, 10, 12, (3D) 7, 11, 13, 14, (4D) 15 (shapes)
    //                 0       1  2  3  4       5  6  7  8   9  10      11  12  13  14       15 (types)
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

    // Set up topstar_offset_ data structure
    topstar_offset_.resize(M);
    for ( Integer i = 0; i < M; ++ i ) {
      topstar_offset_[i] = 0;
      for ( Integer d = 0; d < dimension(); ++ d ) {
        if ( (i & (1L << d)) == 0 ) { 
          topstar_offset_[i] -= place_values_[d];
        }
      }
    }
  }

  /// column
  virtual void
  column ( Integer cell, std::function<void(Integer)> const& callback ) const final {
    Integer shape = cell_shape(cell);
    Integer position = cell % type_size();
    for ( Integer d = 0, bit = 1; d < dimension(); ++ d, bit <<= 1L ) {
      // If cell has no extent in this dimension, no boundaries.
      if ( not (shape & bit) ) continue;
      Integer type_offset = type_size() * ( TS() [ shape ^ bit ] );
      callback( position + type_offset );
      // Otherwise, the cell does have extent in this dimension.
      // It is always the case that such a cell has a boundary to the left.
      Integer right_position = position + PV()[d];
      if (right_position >= type_size()) right_position -= type_size();
      callback( right_position + type_offset );
    }
  }

  /// row
  virtual void
  row ( Integer cell, std::function<void(Integer)> const& callback ) const final {
    Integer shape = cell_shape(cell);
    Integer position = cell % type_size();
    for ( Integer d = 0, bit = 1; d < dimension(); ++ d, bit <<= 1L ) {
      // If cell has extent in this dimension, no coboundaries.
      if ( shape & bit ) continue;
      Integer type_offset = type_size() * ( TS() [ shape ^ bit ] );
      callback( position + type_offset );
      Integer left_position = position - PV()[d];
      if (left_position < 0) left_position += type_size();
      callback( left_position + type_offset );
    }
  }

  /// topstar
  ///   return top dimensional cells in star
  ///   note: assumed twisted periodic conditions
  virtual std::vector<Integer>
  topstar ( Integer cell ) const {
    std::vector<Integer> result;
    Integer shape = cell_shape(cell);
    // Loop through dimension()-bit bitcodes
    Integer M = 1L << dimension(); // i.e. 2^dimension()
    // Compute the topcell x we get by expanding cell to the right in all collapsed dimensions:
    Integer position = cell % type_size(); // (M-1)==(2^D-1) is type of a topcell
    Integer offset = type_size() * (M-1);
    // std::cout << " cell = " << cell << "\n";
    // std::cout << "all-1's topcell = " << x << "\n";
    for ( Integer i = 0; i < M; ++ i ) {
      if ( (shape & ~i ) == 0 ) { // if shape bit is one, then index i bit must be on.
        // i is a valid offset
        // std::cout << " shape = " << shape << " i == " << i << "\n";
        // std::cout << " i = " << i << "\n";
        // std::cout << " topstar_offset_[i] = " << topstar_offset_[i] << "\n";
        result.push_back(offset + (position + topstar_offset_[i] + type_size() ) % type_size());
      } else {
        // TODO skip many invalid indices simultaneously
      }
    }
    return result;
  }

  /// parallelneighbors
  ///   return top dimensional cells in star
  ///   note: assumed twisted periodic conditions
  std::vector<Integer>
  parallelneighbors ( Integer cell ) const {
    Integer shape = cell_shape(cell);
    Integer position = cell % type_size();
    Integer type_offset = type_size() * TS() [ shape ];
    std::vector<Integer> result;
    auto record = [&](Integer offset) {
      Integer k = position + offset;
      if ( k >= type_size()) k -= type_size(); else if ( k < 0 ) k += type_size();
      result.push_back(k + type_offset);
    };
    auto D = dimension();
    std::vector<int> x(D, -1);
    Integer offset = 0;
    for ( Integer d = 0, bit = 1; d < D; ++ d, bit <<= (Integer) 1 ) if ( shape & bit ) offset -= PV()[d];
    record(offset);
    while (1) {
      for ( Integer d = 0, bit = 1; d <= D; ++ d, bit <<= (Integer) 1 ) {
        if ( d == D ) return result;
        if ( shape & bit ) {
          if ( ++x[d] == 2 ) { 
            x[d] = -1; 
            offset -= 2*PV()[d];
          } else { 
            offset += PV()[d];
            record(offset);
            break;
          }
        }
      }
    }
  }

  /// Features

  /// left
  ///   Give cell to "left" in given dimension. 
  ///   Note: uses "twisted" periodic boundary conditions (inconsistent with periodic and acyclic conditions)
  Integer
  left ( Integer cell, Integer dim ) const {
    Integer shape = cell_shape(cell);
    Integer bit = ((Integer)1) << dim;
    Integer position = cell % type_size();
    Integer type_offset = type_size() * TS() [ shape ^ bit ];
    if ( not (shape & bit) ) position -= PV()[dim];
    if (position < 0) position += type_size();
    return position + type_offset;
  }

  /// right
  ///   Give cell to "right" in given dimension. 
  ///   Note: uses "twisted" periodic boundary conditions (inconsistent with periodic and acyclic conditions)
  Integer
  right ( Integer cell, Integer dim ) const {
    Integer shape = cell_shape(cell);
    Integer bit = ((Integer)1) << dim;
    Integer position = cell % type_size();
    Integer type_offset = type_size() * TS() [ shape ^ bit ];
    if ( (shape & bit) ) position += PV()[dim];
    if (position >= type_size()) position -= type_size();
    return position + type_offset;
  }

  /// leftfringe
  ///   cells for which left coboundary wraps around if periodic
  bool
  leftfringe ( Integer cell ) const {
    Integer shape = cell_shape(cell);
    lldiv_t coordinate = {(int64_t)cell, 0}; // (quotient, remainder), see std::div
    for ( Integer d = 0, bit = 1; d < dimension(); ++ d, bit <<= 1L ) {
      coordinate = std::lldiv(static_cast<Integer>(coordinate.quot), boxes()[d] ); 
      if ( (shape & bit) ) continue;
      if ( coordinate.rem == 0 ) return true;
    }
    return false;
  }

  /// rightfringe
  ///    Return true if right boundary would wrap around if periodic 
  bool
  rightfringe ( Integer cell ) const {
    Integer shape = cell_shape(cell);
    lldiv_t coordinate = {(int64_t)cell, 0}; // (quotient, remainder), see std::div
    for ( Integer d = 0, bit = 1; d < dimension(); ++ d, bit <<= 1L ) {
      coordinate = std::lldiv(static_cast<Integer>(coordinate.quot), boxes()[d] ); 
      if ( not (shape & bit) ) continue;
      if ( coordinate.rem + 1 == boxes()[d]) return true;
    }
    return false;
  }

  /// mincoords
  ///   Return integer with 2^i bit set for each dimension i for which 
  ///   cell is at maximum coordinate for dimension i
  Integer
  mincoords ( Integer cell ) const {
    Integer result = 0;
    lldiv_t coordinate = {(int64_t)cell, 0}; // (quotient, remainder), see std::div
    for ( Integer d = 0, bit = 1; d < dimension(); ++ d, bit <<= 1L ) {
      coordinate = std::lldiv(static_cast<Integer>(coordinate.quot), boxes()[d] ); 
      if ( coordinate.rem == 0 ) result |= bit;
    }
    return result;
  }

  /// maxcoords
  ///   Return integer with 2^i bit set for each dimension i for which 
  ///   cell is at maximum coordinate for dimension i
  Integer
  maxcoords ( Integer cell ) const {
    Integer result = 0;
    lldiv_t coordinate = {(int64_t)cell, 0}; // (quotient, remainder), see std::div
    for ( Integer d = 0, bit = 1; d < dimension(); ++ d, bit <<= 1L ) {
      coordinate = std::lldiv(static_cast<Integer>(coordinate.quot), boxes()[d] ); 
      if ( coordinate.rem + 1 == boxes()[d]) result |= bit;
    }
    return result;
  }

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

  /// barycenter
  ///   Give integer barycenter ( doubled coordinates with +1 for directions with extent)
  std::vector<Integer>
  barycenter ( Integer cell ) const {
    auto result = coordinates(cell);
    auto shape = cell_shape(cell);
    for ( Integer i = 0, bit = 1; i < dimension(); ++ i, bit <<= (Integer) 1 ) {
      result[i] <<= (Integer) 1;
      if ( shape & bit ) result[i] += 1;
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

  /// cell_pos
  Integer
  cell_pos ( Integer cell ) const {
    return cell % type_size();
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

  /// cell_dim
  ///   Return dimension of cell
  Integer
  cell_dim ( Integer cell ) const {
    return popcount_(cell_shape(cell));
  }

  /// operator ==
  bool
  operator == ( CubicalComplex const& rhs ) const {
    return boxes() == rhs.boxes();
  }

  /// operator <
  bool
  operator < ( CubicalComplex const& rhs ) const {
    return std::lexicographical_compare(boxes().begin(), boxes().end(),
                                          rhs.boxes().begin(), rhs.boxes().end());
  }

  /// operator <<
  friend std::ostream & operator << ( std::ostream & stream, CubicalComplex const& stream_me ) {
    stream << "CubicalComplex([";
    for ( auto x : stream_me.boxes() ) stream << x << ",";
    return stream << "])";
  }

  /// print_cell (for debugging )
  void
  print_cell (Integer cell_index) const {
    std::cout << cell_index << "\n";
    std::cout << "  coordinates(" << cell_index << ") = [";
    for ( auto x : coordinates(cell_index) ) std::cout << x << ", "; std::cout << "]\n";
    std::cout << "  shape(" << cell_index << ") = " << cell_shape(cell_index) << "\n";
  }

  Integer
  type_size ( void ) const {
    return type_size_;
  }

  /// TS
  ///   Given shape, return type
  std::vector<Integer> const&
  TS ( void ) const {
    return type_from_shape_;
  }

  /// TS
  ///   Given type, return shape
  std::vector<Integer> const&
  ST ( void ) const {
    return shape_from_type_;
  }

  /// PV
  ///   Return 
  std::vector<Integer> const&
  PV ( void ) const {
    return place_values_;
  }

private:

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
  std::vector<Integer> topstar_offset_;
  Integer num_types_;
  Integer type_size_;
};

/// std::hash<CubicalComplex>
namespace std {
  template<> struct hash<CubicalComplex> {
    typedef CubicalComplex argument_type;
    typedef std::size_t result_type;
    result_type operator()(argument_type const& complex) const {
      using pychomp::hash_value;
      using pychomp::hash_combine;
      std::size_t seed = 0;
      for ( auto x : complex.boxes() ) {
        hash_combine(seed,hash_value(x));
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
    .def("boxes", &CubicalComplex::boxes)
    .def("coordinates", &CubicalComplex::coordinates)
    .def("barycenter", &CubicalComplex::barycenter)    
    .def("cell_type", &CubicalComplex::cell_type)
    .def("cell_shape", &CubicalComplex::cell_shape)
    .def("cell_pos", &CubicalComplex::cell_pos)
    .def("cell_dim", &CubicalComplex::cell_dim)
    .def("cell_index", &CubicalComplex::cell_index)
    .def("left", &CubicalComplex::left)
    .def("right", &CubicalComplex::right)
    .def("leftfringe", &CubicalComplex::leftfringe)
    .def("rightfringe", &CubicalComplex::rightfringe)
    .def("mincoords", &CubicalComplex::mincoords)
    .def("maxcoords", &CubicalComplex::maxcoords)
    .def("parallelneighbors", &CubicalComplex::parallelneighbors);
}
