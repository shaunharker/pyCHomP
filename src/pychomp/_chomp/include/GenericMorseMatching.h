/// GenericMorseMatching.h
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

class GenericMorseMatching : public MorseMatching {
public:
  /// GenericMorseMatching
  GenericMorseMatching ( std::shared_ptr<Complex> complex_ptr ) {
    std::shared_ptr<GradedComplex> graded_complex_ptr (new GradedComplex(complex_ptr, [](Integer i){return 0;}));
    construct(graded_complex_ptr);
  }

  /// GenericMorseMatching
  // DRY mistake -- only a few lines differ. 
  GenericMorseMatching ( std::shared_ptr<GradedComplex> graded_complex_ptr ) {
    construct(graded_complex_ptr);
  }

  /// construct
  void
  construct ( std::shared_ptr<GradedComplex> graded_complex_ptr ) {
    GradedComplex const& graded_complex = *graded_complex_ptr;
    Complex const& complex = *graded_complex.complex();
    Integer N = complex.size();
    mate_.resize(N,-1);
    priority_.resize(N);
    Integer num_processed = 0;
    std::vector<Integer> boundary_count (N);
    std::unordered_set<Integer> coreducible;
    std::unordered_set<Integer> ace_candidates;

    auto bd = [&](Integer x) {
      Chain result;
      auto x_val = graded_complex.value(x);
      for ( auto y : complex.boundary({x}) ) {
        auto y_val = graded_complex.value(y);
        if ( y_val > x_val ) {
          throw std::logic_error("graded_complex closure property failed line MorseMatching line 98");
        }
        if ( x_val == y_val ) result += y;
      }
      return result;
    };

    auto cbd = [&](Integer x) {
      Chain result;
      auto x_val = graded_complex.value(x);
      for ( auto y : complex.coboundary({x}) ) {
        if ( x_val == graded_complex.value(y) ) result += y;
      }
      return result;
    };

    for ( auto x : complex ) {
      boundary_count[x] = bd(x).size();
      switch ( boundary_count[x] ) {
        case 0: ace_candidates.insert(x); break;
        case 1: coreducible.insert(x); break;
      }
    }

    auto process = [&](Integer y){
      priority_[y] = graded_complex.value(y)*complex.size() + num_processed ++;
      coreducible.erase(y);
      ace_candidates.erase(y);
      for ( auto x : cbd(y) ) {
        boundary_count[x] -= 1;
        switch ( boundary_count[x] ) {
          case 0: coreducible.erase(x); ace_candidates.insert(x); break;
          case 1: coreducible.insert(x); break;
        }
      }
    };

    while ( num_processed < N ) {
      if ( not coreducible.empty() ) {
        Integer K, Q;
        // Extract K
        auto it = coreducible.begin(); K = *it; coreducible.erase(it); // pop from unordered_set
        // Find mate Q
        for ( auto x : bd(K) ) if ( mate_[x] == -1 ) { Q = x; break; }
        if ( graded_complex.value(K) != graded_complex.value(Q) ) {
          throw std::logic_error("graded_complex error -- memory bug? MorseMatching line 132");
        }
        mate_[K] = Q; mate_[Q] = K;
        process(Q); process(K);
      } else {
        Integer A;
        // Error: what if there are zero ace candidates?
        auto it = ace_candidates.begin(); A = *it; ace_candidates.erase(it); // pop from unordered_set
        mate_[A] = A;
        process(A);
      }
    }

    // Compute critical cells
    Integer D = complex.dimension();
    //std::cout << "  Dimension = " << dim_ << "\n";
    begin_.resize(D+2);
    Integer idx = 0;
    for ( Integer d = 0; d <= D; ++ d ) {
      begin_[d] = idx;
      //std::cout << "  begin_[" << d << "] = " << *begin_[d] << "\n";
      for ( auto v : complex(d) ) {
        //std::cout << "  Inspecting cell " << v << " with mate " << matching_ -> mate(v) << "\n";
        if ( mate(v) == v ) { 
          //std::cout << "  Identified cell " << idx << "\n";
          reindex_.push_back({v, idx});
          ++idx;
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
    return mate_[x];
  }

  /// priority
  Integer
  priority ( Integer x ) const { 
    return priority_[x];
  }

private:
  std::vector<Integer> mate_;
  std::vector<Integer> priority_;
  BeginType begin_;
  ReindexType reindex_;
};

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline void
GenericMorseMatchingBinding(py::module &m) {
  py::class_<GenericMorseMatching, std::shared_ptr<GenericMorseMatching>>(m, "GenericMorseMatching")
    .def(py::init<std::shared_ptr<Complex>>())
    .def(py::init<std::shared_ptr<GradedComplex>>())    
    .def("mate", &GenericMorseMatching::mate)
    .def("priority", &GenericMorseMatching::priority);
}
