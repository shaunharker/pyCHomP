/// CubicalComplex.h
/// Shaun Harker 2016-12-15-2159
/// MIT LICENSE

#pragma once

#include "common.h"

#include "Integer.h"
#include "Iterator.h"
#include "Chain.h"
#include "Complex.h"

///   CubicalComplex  TwistedCubicalComplex     RightWrap
///                      | # | #                | # | #
///      * - *     i     * - * -        p             -    
/// 0--> | # |   c--->   | # | #      ---->>          #   ---> 0      
///      * - *           * - * -                      -         
///
///   This is a short exact sequence. 
///

///  Rank-Select 
///    Given a sequence A of 0's and 1's of length n, 
///    i.e. A \in {0,1}^n, we consider two types of queries:
///      1. rank query:  
///         rank : {0,1}^n x [n] -> [n]
///         rank(A, i) := 
///           #{ k \in [n] : k < i and A_k = 1}
///      2. select query:
///         select : {0,1}^n x [n] -> [n+1]
///         select(A, i) := 
///           \max { k \in [n+1] : rank(A, k) <= i}
///
///    012345678  indices
///    000100010  data
///    000111122  rank
///    379999999  select
/// 
///    f and g are Galois connection =
///      f x <= y  iff x <= g y  
///      
///   claim (rank(A,-), select(A,-)) is a galois connection
///      x 012345678       y 012345678
///   rank 000111122  select 379999999
///    
///   select(A, 1) == 7
///   rank(A, 7)   == 2
///   select(A, 1) == 
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
    // needs to give wrap to boxes
    boxes_ = boxes;
    std::vector<Integer> wrapped_boxes = boxes;
    for ( auto & x : boxes ) ++ x;
    twisted_complex.assign(wrapped_boxes);
    //  Integer dim_;
    //  std::vector<Iterator> begin_;
    // needs to set these
    dim_ = boxes.size();
    begin_.resize(dim_+2);
    for ( int d = 0; d < dim_+1; ++ d) {
      begin_[d] = Iterator(tidx_to_cidx(*twisted_complex.begin(d)));
    }
    begin_[dim_+1] = // total number of cells;

  }

  /// boundary
  virtual Chain
  boundary ( Chain const& chain ) const {
    return project(twisted_complex.boundary(include(chain)));
  }

  /// coboundary
  virtual Chain
  coboundary ( Chain const& chain ) const {
    return project(twisted_complex.coboundary(include(chain)));
  }

  /// topstar
  ///   return top dimensional cells in star
  ///   note: assumed twisted periodic conditions
  virtual std::vector<Integer> // should this use Chain?
  topstar ( Integer cell ) const {
    return project_vec(twisted_complex.topstar(cidx_to_tidx(cell)));
  }

  /// parallelneighbors
  ///  Gives cells in star of closure with same shape
  std::vector<Integer>
  parallelneighbors ( Integer cell ) const {
    return project_vec(twisted_complex.parallelneighbors(cidx_to_tidx(cell)));
  }

  /// Features

  /// left
  ///   Give cell to "left" in given dimension. 
  ///   Return -1 if no such cell exists
  Integer
  left ( Integer cell, Integer dim ) const {
    return tidx_to_cidx(twisted_complex.left(cidx_to_tidx(cell)));
  }

  /// right
  ///   Give cell to "right" in given dimension. 
  ///   Return -1 if no such cell exists
  Integer
  right ( Integer cell, Integer dim ) const {
    return tidx_to_cidx(twisted_complex.right(cidx_to_tidx(cell)));
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
    return twisted_complex.coordinates(cidx_to_tidx(cell));
  }

  /// barycenter
  ///   Give integer barycenter ( doubled coordinates with +1 for directions with extent)
  std::vector<Integer>
  barycenter ( Integer cell ) const {
    return twisted_complex.barycenter(cidx_to_tidx(cell));
  }

  Integer
  cell_type ( Integer cell ) {
    return twisted_complex.cell_type(cidx_to_tidx(cell));
  }

  /// cell_shape
  ///  Returns Integer with ith bit set if cell
  ///  has non-zero width in dimension i
  ///   i.e. = \sum_{i=0}^D 2^i * X(i)
  ///   where X(i) = 1 if (cell has with in dim i) else 0
  Integer
  cell_shape ( Integer cell ) {
    return twisted_complex.cell_shape(cidx_to_tidx(cell));
  }

  // /// cell_pos
  // Integer
  // cell_pos ( Integer cell ) const {
  //   return twisted_complex.cell_pos(cidx_to_tidx(cell));
  // }

  /// cell_index
  Integer
  cell_index ( std::vector<Integer> const& coordinates, 
               Integer shape ) {
    return tidx_to_cidx(twisted_complex.cell_index(coordinates, shape));
  }

  /// cell_dim
  ///   Return dimension of cell
  Integer
  cell_dim ( Integer cell ) const {
    return twisted_complex.cell_dim(cidx_to_tidx(cell));
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

private:
  std::vector<Integer> boxes_;

  auto 
  cidx_to_tidx
    ( Integer cidx ) 
    ->
    Integer
  { 
    // TODO : create this formula
  }

  auto 
  tidx_to_cidx 
    ( Integer tidx ) 
    ->
    Integer
  { 
    // TODO : create this formula
    // cidx = ...
    if (cidx_to_tidx(cidx) == tidx) {
      return cidx;
    } else {
      return -1;
    }
  }

  auto 
  include(Chain cc) -> Chain
  {
    Chain result;
    for ( auto cidx : cc ) {
      result += {cidx_to_tidx(cidx)};
    }
    return result;
  }

  auto
  valid_twisted_cell(Integer tidx)
    -> Integer 
  {
    return cidx_to_tidx(tidx_to_cidx(tidx)) == tidx;
  }

  auto 
  project_vec(std::vector<Integer> v) 
    -> std::vector<Integer>
  {
    std::vector<Integer> result;
    for ( auto tidx : v ) { 
      auto cidx = tidx_to_cidx(tidx);
      if ( cidx != -1 ) result.push_back(cidx);
    }
    return result;
  }

  auto 
  project(Chain tc) -> Chain
  {
    std::vector<Integer> v ( tc.begin(), tc.end() );
    auto proj_v = project_vec(v);
    return Chain(proj_v.begin(), proj_v.end());
  }
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
    //.def("cell_pos", &CubicalComplex::cell_pos)
    .def("cell_dim", &CubicalComplex::cell_dim)
    .def("cell_index", &CubicalComplex::cell_index)
    .def("left", &CubicalComplex::left)
    .def("right", &CubicalComplex::right)
    .def("parallelneighbors", &CubicalComplex::parallelneighbors);
}
