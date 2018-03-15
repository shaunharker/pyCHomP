/// SimplicialComplex.h
/// Shaun Harker
/// MIT LICENSE
/// 2018-03-09

#pragma once

#include "common.h"

typedef std::vector<Integer> Simplex;

inline std::vector<Simplex>
simplex_boundary(Simplex const& s) {
  //std::cout << "bd ["; for ( auto v : s ) std::cout << v << ", "; std::cout << "] = \n";
  std::vector<Simplex> result;
  if ( s.size() > 1 ) {
    for ( Integer i = 0; i < s.size(); ++ i ) {
      Simplex t = s;
      t.erase(t.begin() + i);
      result.push_back(t);
      //std::cout << "  ["; for ( auto v : t ) std::cout << v << ", "; std::cout << "] + \n";
    }
  }
  return result;
}

/// SimplicialComplex
class SimplicialComplex : public Complex {
public:

  /// SimplicialComplex
  SimplicialComplex ( std::vector<Simplex> const& maximal_simplices );

  /// column
  ///   Apply "callback" method to every element in ith column of
  ///   boundary matrix
  virtual void
  column ( Integer i, std::function<void(Integer)> const& callback) const final;

  /// row
  ///   Apply "callback" method to every element in ith row of
  ///   boundary matrix
  virtual void
  row ( Integer i, std::function<void(Integer)> const& callback) const final;

  /// simplex
  ///   Given a cell index, return the associated Simplex
  Simplex
  simplex ( Integer i ) const;

  /// idx
  ///   Given a simplex object, return the associated cell index.
  ///   If simplex not in complex, return -1.
  Integer
  idx ( Simplex const& s ) const;

private:
  std::unordered_map<Simplex, Integer, pychomp::hash<Simplex>> idx_;
  std::vector<Simplex> simplices_;
  std::vector<Chain> bd_;
  std::vector<Chain> cbd_;
  
  /// add_simplex
  bool
  add_simplex ( Simplex const& s );

  /// add_closed_simplex
  void
  add_closed_simplex ( Simplex const& s );
  
};

inline SimplicialComplex::
SimplicialComplex (std::vector<Simplex> const& max_simplices) {
  //std::cout << "SimplicialComplex " << max_simplices.size() << "\n";
  for ( auto s : max_simplices ) add_closed_simplex ( s );
  Integer N = simplices_.size();
  std::sort(simplices_.begin(), simplices_.end(), 
    []( Simplex const& lhs, Simplex const& rhs ){ 
      if (lhs.size() < rhs.size()) return true; 
      if (lhs.size() > rhs.size()) return false;
      return std::lexicographical_compare(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
    });
  idx_.clear();
  for ( Integer i = 0; i < N; ++ i ) idx_[simplices_[i]] = i;
  dim_ = -1;
  bd_.resize(N);
  cbd_.resize(N);
  for ( Integer i = 0; i < N; ++ i ) {
    Simplex const& s = simplices_[i];
    Integer simplex_dim = s.size() - 1;
    // std::cout << "i = " << i << "\n";
    // std::cout << " s =  ["; for ( auto v : s ) std::cout << v << ", "; std::cout << "] + \n";
    // std::cout << " simplex_dim == " << simplex_dim << "\n";
    // std::cout << " dim_ == " << dim_ << "\n";
    if ( simplex_dim > dim_ ) {
      // std::cout << "New dimension\n";
      ++ dim_;
      begin_.push_back(Iterator(i));
      // std::cout << "Pushed " << i << " onto begin_\n";
    }
    Chain c;
    for ( Simplex const& t : simplex_boundary(s) ) c += idx_[t];
    bd_[i] = c;
    //std::cout << "boundary of " << i << " is equal to " << c << "\n";
  }
  begin_.push_back(Iterator(N));
  // std::cout << "Pushed " << N << " onto begin_\n";

  for ( Integer i = 0; i < N; ++ i ) {
    Chain bd = bd_[i];
    for ( Integer j : bd ) {
      cbd_[j] += i;
    }
  }
}

inline Simplex SimplicialComplex::
simplex ( Integer i ) const{
  return simplices_[i];
}

inline Integer SimplicialComplex::
idx ( Simplex const& s ) const {
  auto it = idx_.find(s);
  if ( it == idx_.end() ) return -1;
  return it -> second;
}

inline bool SimplicialComplex::
add_simplex (Simplex const& s) {
  if ( idx(s) == -1 ) {
    idx_[s] = simplices_.size();
    simplices_.push_back(s);
    return true;
  } else {
    return false;
  }
}

inline void SimplicialComplex::
add_closed_simplex (Simplex const& s) {
  std::stack < Simplex > work_stack;
  work_stack.push(s);
  while ( not work_stack.empty() ) {
    auto t = work_stack.top();
    work_stack.pop();
    bool inserted = add_simplex(t);
    if ( inserted ) {
      for ( auto u : simplex_boundary(t) ) {
        work_stack.push(u);
      }
    }
  }
}

inline void SimplicialComplex::
column ( Integer i, std::function<void(Integer)> const& callback ) const { 
  for ( auto x : bd_[i] ) callback(x);
}

inline void SimplicialComplex::
row ( Integer i, std::function<void(Integer)> const& callback ) const {
  for ( auto x : cbd_[i] ) callback(x);
}

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline void
SimplicialComplexBinding(py::module &m) {
  py::class_<SimplicialComplex, std::shared_ptr<SimplicialComplex>, Complex>(m, "SimplicialComplex")
    .def(py::init<std::vector<Simplex> const&>())
    .def("simplex", &SimplicialComplex::simplex)
    .def("idx", &SimplicialComplex::idx);
}
