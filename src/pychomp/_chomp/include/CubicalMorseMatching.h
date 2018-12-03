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
#include "GradedComplex.h"
#include "MorseMatching.h"

class CubicalMorseMatching : public MorseMatching {
public:
  /// CubicalMorseMatching
  CubicalMorseMatching ( std::shared_ptr<CubicalComplex> complex_ptr ) : complex_(complex_ptr) {
    type_size_ = complex_ -> type_size();
    graded_complex_.reset(new GradedComplex(complex_, [](Integer i){return 0;}));
  }

  /// idea 1. Pass an an argument to the constructor of CubicalMorseMatching
  ///         that abstracts the parallelism capability
  /// idea 2. Make a new constructor that glues together results for partial objects
  ///  conj. the induced matching on the subdivided graded complexes are the subdivided induced matchings
  ///       * we can subdivide in a way that respects the implicit matching
  ///   * A "left-closed" subcomplex  *-----*------
  ///   * how to define the grading on left-closed subcomplexes?
  ///   (1) for a complex with values on top cells, we can induce a grading.
  ///   problem. If we take a left-closed subcomplex, we will induce a different grading using the top cells
  ///            in the subcomplex.
  ///   So we need extra information so the induced grading on a LC-subcomplex
  ///   Another problem: the twisted periodicity implementation could result in strange behaviors.
  ///    Two birds with one stone: can we fix the twisted periodicity problem with fringe and put values on that
  ///    fringe to make the values on the left-fringe correct?
  /// (A) top-cell valuations induce cell gradations
  /// Types:
  ///   CubicalComplex  --- this is really a cubical complex with twisted periodicity.
  ///   CubicalMorseMatching
  ///   GradedComplex
  ///   We don't have: GradedCubicalComplex
  ///                  GradedCubicalComplexWhereGradingIsInducedFromTopCellValuation
  ///                  CubicalComplexWithTopCellValuation
  ///                  CubicalComplexWithTopCellValuationAndFringeData
  ///   we have: construct a grading from a top-cell valuation on complex
  ///   we need: construct a grading on a left-closed subcomplex given a top-cell valuation on complex
  ///            idea. implement the left-closed subcomplex with a cubicalcomplex with a wrapper
  ///                  compute a top-cell valuation on the underlying-implementation-complex
  ///                  compute a grading on the underlying-implementation-complex with this top-cell valuation
  ///                  give a view to the 
  /// ...
  ///   Start with a CubicalComplex with a top cell valuation
  ///   Take a subcomplex out of it with a left-wrapper, and also a sub top-cell valuation correspond to it
  ///   Now induce a grading on the subcomplex.
  ///   Now take a sub-gradedcomplex by removing left-wrapper.
  /// Plan. 1. Refactor CubicalComplex so that twisted periodicity is an implementation detail
  ///          (Say, by making a TwistedCubicalComplex class and making a CubicalComplex which uses it
  ///           an as underlying implementation -- i.e. a subcomplex of Twisted with right-fringe missing)
  ///       2. (TwistedCubicalComplex, TopCellValuation) -> GradedComplex
  ///       3. class LeftExtendedTopCellValuation
  ///       4. class (CubicalComplex, LeftExtendedTopCellValuation)
  ///          -- implemented with (TwistedCubicalComplex, TopCellValuation)
  ///       5. subdivision procedure for (CubicalComplex, LeftExtendedTopCellValuation)
  ///
  ///   Alternative design with grading as fundamental:
  ///       1. (TwistedCubicalComplex, GradingFromTopCellValuation)
  ///       2. (CubicalComplex, GradingInducedFromLeftExtendedTopCellValuation) 
  ///            -- impl as (TwistedCubicalComplex, GradingFromTopCellValuation)
  ///       3. subdivide : 
  ///            (CubicalComplex, GradingInducedFromLeftExtendedTopCellValuation) x dim -> 
  ///              (CubicalComplex, GradingInducedFromLeftExtendedTopCellValuation) x 
  ///              (CubicalComplex, GradingInducedFromLeftExtendedTopCellValuation)
  ///       4. subdivide1 :  CubicalComplex x dim -> CubicalComplex x CubicalComplex
  ///       5. subdivide2 :  GIFLETCV x dim -> GIFLETCV x GIFLETCV
  ///       6. subdivide (x , y) dim = (subdivide1 x dim , subdivide2 y dim)
  ///       7. glue: CubicalComplex x CubicalComplex x dim -> CubicalComplex
  ///          glue: GIFLETCV x GIFLETCV x dim -> GIFLETCV
  ///     actually, for subdivide, we don't want to just get a pair of cubicalcomplex.
  ///     i d = d i won't hold.   * ------ *
  ///                             * ------
  ///     p d = d p ? 
  /// 
  ///      * -- *     i     * - *    p           
  /// 0-->      |   c--->   | # |  ---->> | #  ---> 0      
  ///           *           * - *         * -          
  ///
  ///
  ///     if subcomplex means a subobject in the category of chain complexes,
  ///     then the formal definition is (S, i : S to C) 
  ///     (Q, p : C to Q)
  ///
  ///         p          
  ///     C ---->> Q    every epi in an abelian category is a cokernel.
  ///                  
  ///     left-closed quotient complexes.
  ///     axioms for subdivision: 
  ///
  ///
  ///           |#          | # | #         | #
  ///           *-    i     * - * -    p    * -       
  /// 0-->      |#  c--->   | # | #  ---->> | #  ---> 0      
  ///           *-          * - * -         * -          
  ///
  ///   subdivide : CubicalComplex -> QuotientComplex x SubComplex
  ///   glue      : QuotientComplex x SubComplex -> CubicalComplex
  ///   

  /// CubicalMorseMatching
  CubicalMorseMatching ( std::shared_ptr<GradedComplex> graded_complex_ptr ) : graded_complex_(graded_complex_ptr) {
    complex_ = std::dynamic_pointer_cast<CubicalComplex>(graded_complex_->complex());
    if ( not complex_ ) {
      throw std::invalid_argument("CubicalMorseMatching must be constructed with a Cubical Complex");
    }
    type_size_ = complex_ -> type_size();
    Integer D = complex_ -> dimension();
    Integer idx = 0;
    begin_.resize(D+2);
    for ( Integer d = 0; d <= D; ++ d) {
      begin_[d] = idx;
      for ( auto v : (*complex_)(d) ) { // TODO: skip fringe cells
        if ( ! complex_ -> rightfringe(v) ) {
          if ( mate(v) == v ) { 
            reindex_.push_back({v,idx});
            ++idx;
          }
        }
      }
    }
    begin_[D+1] = idx;
  }

  /// critical_cells
  std::pair<BeginType const&,ReindexType const&>
  critical_cells ( void ) const {
    return {begin_,reindex_};
  }

  /// mate
  Integer
  mate ( Integer x ) const { 
    return mate_(x, complex_ -> dimension());
  }

  /// priority
  Integer
  priority ( Integer x ) const { 
    return type_size_ - x % type_size_;
  }

private:
  uint64_t type_size_;
  std::shared_ptr<GradedComplex> graded_complex_;
  std::shared_ptr<CubicalComplex> complex_;
  BeginType begin_;
  ReindexType reindex_;

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
  // Note: the reason for the "fringe" check preventing mating is that otherwise it is possible to 
  //       end up with a cycle 
  // TODO: Furnish a proof of correctness and complexity that this cannot produce cycles.
  Integer mate_ ( Integer cell, Integer D ) const {
    //bool fringe = complex_ -> rightfringe(cell);
    if ( complex_ -> rightfringe(cell) ) return cell; // MAYBE
    //Integer mincoords = complex_ -> mincoords(cell); // TODO: optimize to compute this as it loops through d rather than demanding all
    //Integer maxcoords = complex_ -> maxcoords(cell); // TODO: optimize to compute this as it loops through d rather than demanding all
    Integer shape = complex_ -> cell_shape(cell);
    Integer position = cell % complex_ -> type_size();
    if ( position == complex_ -> type_size() - 1 ) return cell; // Break cycles
    for ( Integer d = 0, bit = 1; d < D; ++ d, bit <<= 1L  ) {
      // If on right fringe for last dimension, prevent mating with left boundary
      if ( (d == D-1) && (position + complex_ -> PV()[d] >= complex_ -> type_size()) ) break;
      //if ( fringe && (mincoords & bit) ) continue; // Todo -- is this the best
      //if ( bit & maxcoords ) continue; // Don't connect fringe to acyclic part
      Integer type_offset = complex_ -> type_size() * complex_ -> TS() [ shape ^ bit ];
      Integer proposed_mate = position + type_offset;
      if ( graded_complex_ -> value(proposed_mate) == graded_complex_ -> value(cell) && proposed_mate == mate_(proposed_mate, d) ) { 
        return proposed_mate;
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
    .def(py::init<std::shared_ptr<GradedComplex>>())    
    .def("mate", &CubicalMorseMatching::mate)
    .def("priority", &CubicalMorseMatching::priority);
}
