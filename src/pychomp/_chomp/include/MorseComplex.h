/// MorseComplex.h
/// Shaun Harker
/// 2017-07-19
/// MIT LICENSE

#pragma once

#include <memory>
#include <queue>
#include <tuple>
#include <unordered_set>
#include <vector>

#include "Integer.h"
#include "Iterator.h"
#include "Chain.h"
#include "Complex.h"
#include "MorseMatching.h"

class MorseComplex : public Complex {
public:


  /// MorseComplex constructor
  MorseComplex ( std::shared_ptr<Complex> arg_base, 
                 std::shared_ptr<MorseMatching> arg_matching ) 
               : base_(arg_base), matching_(arg_matching) {
    //std::cout << "MorseComplex constructor\n";
    dim_ = base()->dimension();
    //std::cout << "  Dimension = " << dim_ << "\n";
    begin_.resize(dim_+2);
    Integer idx = 0;
    for ( Integer d = 0; d <= dim_; ++ d ) {
      begin_[d] = Iterator(idx);
      //std::cout << "  begin_[" << d << "] = " << *begin_[d] << "\n";
      for ( auto v : (*base())(d) ) {
        //std::cout << "  Inspecting cell " << v << " with mate " << matching_ -> mate(v) << "\n";
        if ( matching_ -> mate(v) == v ) { 
          //std::cout << "  Identified cell " << idx << "\n";
          project_[v] = idx++;
          include_.push_back(v);
        }
      }
    }
    begin_[dim_+1] = idx;
    //std::cout << "  Total number of cells: " << idx << "\n";

    // boundary
    bd_.resize(size());
    //std::cout << "MorseComplex. There are " << size() << " cells.\n";
    for ( auto ace : *this ) {
      //std::cout << "  Computing boundary for cell ace ==" << ace << "\n";
      //std::cout << "     include({ace}) = " << include({ace}) << "\n";
      bd_[ace] = lower(base()->boundary(include({ace})));
      //std::cout << "     bd(ace) = " << bd_[ace] << "\n";
    }

    // coboundary
    cbd_.resize(size());
    for ( auto ace : *this ) {
      for ( auto bd_cell : bd_[ace] ) cbd_[bd_cell] += ace;
    }

  }

  /// delegating constructor
  MorseComplex ( std::shared_ptr<Complex> arg_base ) 
               : MorseComplex(arg_base, MorseMatching::compute_matching(arg_base)) {

  }

  /// boundary
  virtual Chain
  boundary ( Chain const& c ) const final {
    Chain result;
    for ( auto x : c ) result += bd_[x];
    return result;
  }

  /// coboundary
  virtual Chain
  coboundary ( Chain const& c ) const final {
    Chain result;
    for ( auto x : c ) result += cbd_[x];
    return result;
  }

  /// column
  ///   Apply "callback" method to every element in ith column of
  ///   boundary matrix
  virtual void
  column ( Integer i, std::function<void(Integer)> const& callback) const final {
    for ( auto x : bd_[i] ) callback(x);
  };

  /// row
  ///   Apply "callback" method to every element in ith row of
  ///   boundary matrix
  virtual void
  row ( Integer i, std::function<void(Integer)> const& callback) const final {
    for ( auto x : cbd_[i] ) callback(x);
  };
  

  // Feature

  /// base
  std::shared_ptr<Complex>
  base ( void ) const {
    return base_;
  }

  /// matching
  std::shared_ptr<MorseMatching>
  matching ( void ) const {
    return matching_;
  }

  /// include
  Chain
  include ( Chain const& c ) {
    Chain result;
    for ( auto x : c ) result += include_[x];
    return result;
  }

  /// project
  Chain
  project ( Chain const& c ) {
    Chain result;
    for ( auto x : c ) { 
      if ( project_.count(x) > 0 ) result += project_[x];
    }
    return result;
  }

  /// lift
  Chain
  lift ( Chain const& c ) {
    Chain included = include ( c );
    Chain canonical; Chain gamma;
    std::tie(canonical, gamma) = flow ( base () -> boundary ( included ) );
    return included + gamma;
  }

  /// lower
  Chain
  lower ( Chain const& c ) {
    Chain canonical; Chain gamma;
    std::tie(canonical, gamma) = flow ( c );
    return project(canonical);
  }

  /// flow
  std::pair<Chain, Chain>
  flow ( Chain const& input ) const {
    //std::cout << "MorseComplex::flow\n";
    Chain canonical, gamma;
    std::unordered_set<Integer> queens;
    auto compare = [&](Integer x, Integer y){return matching_ -> priority(x) < matching_ -> priority(y);};
    std::priority_queue<Integer, std::vector<Integer>, decltype(compare)> priority ( compare );
    auto isQueen = [&](Integer x){ return x < matching_ -> mate(x); };

    // auto process = [&](Chain const& c ) {
    //   for ( auto x : c ) {
    //     if ( isQueen(x) && queens.count(x) == 0) {
    //       queens . insert (x);
    //       priority . push (x);
    //     }
    //   }
    //   canonical += c;
    // };

    auto process = [&](Integer x) {
      if ( isQueen(x) && queens.count(x) == 0) {
        queens . insert (x);
        priority . push (x);
      }
      canonical += x;
    };

    for ( auto x : input ) process(x);

    while ( not priority . empty () ) {
      //std::cout << "  Current chain = " << canonical << "\n";
      auto queen = priority.top(); priority.pop();
      if ( canonical . count ( queen ) == 0 ) continue;
      auto king = matching_ -> mate ( queen );
      gamma += king;
      //std::cout << "    Reducing queen " << queen << " with king " << king << " and priority " << matching_->priority(queen) << "\n";
      //std::cout << "       The boundary of king is " << base()->boundary({king})<<"\n";
      base() -> column(king, process);
      //process( base()->boundary({king}) );
    }
    //std::cout << "  COMPLETE chain = " << canonical << "\n";

    return std::tie(canonical, gamma); // TODO -- prevent copy? or optimizer already does?
  }

private:
  std::shared_ptr<Complex> base_;
  std::shared_ptr<MorseMatching> matching_;
  std::vector<Integer> include_;
  std::unordered_map<Integer, Integer> project_;
  std::vector<Chain> bd_;
  std::vector<Chain> cbd_;
};


/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline void
MorseComplexBinding(py::module &m) {
  py::class_<MorseComplex, std::shared_ptr<MorseComplex>, Complex>(m, "MorseComplex")
    .def(py::init<std::shared_ptr<Complex>, std::shared_ptr<MorseMatching>>())
    .def(py::init<std::shared_ptr<Complex>>())
    .def("include", &MorseComplex::include)
    .def("project", &MorseComplex::project)
    .def("lift", &MorseComplex::lift)
    .def("lower", &MorseComplex::lower)
    .def("flow", &MorseComplex::flow)
    .def("base", &MorseComplex::base)
    .def("matching", &MorseComplex::matching);
}
