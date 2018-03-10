/// SimplicialComplex.h
/// Shaun Harker
/// MIT LICENSE
/// 2018-03-09

#pragma once

typedef std::vector<Integer> Simplex;

inline
std::size_t hash_value(Simplex const& simplex ) {  
  std::size_t seed = 0;
  for ( auto v : simplex ) seed += pychomp::hash_combine ( seed, simplex [ i ] );
  return seed;
}

inline std::vector<Simplex>
simplex_boundary(Simplex const& s) {
  std::vector<Simplex> result;
  if ( s.size() > 1 ) {
    for ( Integer i = 0; i < s.size(); ++ i ) {
      Simplex t = s;
      t.erase(t.begin() + i);
      result.push_back(t);
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
  column ( Integer i, std::function<void(Integer)> const& callback) const;

  /// row
  ///   Apply "callback" method to every element in ith row of
  ///   boundary matrix
  virtual void
  row ( Integer i, std::function<void(Integer)> const& callback) const;

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
  std::unordered_map<Simplex, Integer> idx_;
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


inline void SimplicialComplex::
SimplicialComplex (std::vector<Simplex> const& max_simplices) {
  for ( auto s : max_simplices ) add_closed_simplex ( s );
  Integer N = simplices_.size();
  std::sort(simplices_.begin(), simplices_.end(), []( Simplex const& lhs, Simplex const& rhs ){ return lhs.size() < rhs.size(); });
  idx_.clear();
  for ( Integer i = 0; i < N; ++ i ) idx_[simplices_[i]] = i;
  dim_ = -1;
  bd_.resize(N);
  cbd_.resize(N);
  for ( Integer i = 0; i < N; ++ i ) {
    Simplex const& s = simplices_[i];
    auto simplex_dim = s.size() - 1;
    if ( simplex_dim > dim_ ) {
      ++ dim_;
      begin_.push_back(Iterator(i));
    }
    Chain c;
    for ( Simplex const& t : simplex_boundary(s) ) c += idx_[t];
    bd_[i] = c;
  }
  begin_.push_back(Iterator(simplices_.size()));


  // Fill in coboundary data
  for ( Index i = 0; i < size ( dim ); ++ i ) {
    Chain bd = boundary ( i, dim );
    for ( Term const& t : bd () ) {
      coboundaries_ [ dim - 1 ] [ t . index () ] += Term ( i,  t . coef () );
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
  if ( s == idx_.end() ) return -1;
  return *it;
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
  for ( auto x : bd_[i] ) callback(i);
}

inline void SimplicialComplex::
row ( Integer i, std::function<void(Integer)> const& callback ) const {
  for ( auto x : cbd_[i] ) callback(i);
}
