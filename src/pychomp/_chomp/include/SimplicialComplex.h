/// SimplicialComplex.h
/// Shaun Harker
/// MIT LICENSE
/// 2018-03-09

#pragma once

class Simplex {
private:
  std::vector<Integer> simplex_;
  
public:
  /// Simplex
  Simplex ( void ) {};

  /// Simplex
  Simplex (std::vector<Integer> const& simplex ) : simplex_ ( simplex ) {}

  /// simplex
  std::vector < Integer > const&
  simplex ( void ) const { return simplex_; }

  /// simplex
  std::vector < Integer > &
  simplex ( void ) { return simplex_; }
  
  /// dimension
  unsigned Integer 
  dimension ( void ) const { return simplex_ . size () - 1; }

  /// operator []
  Integer
  operator [] ( Integer i ) const { return simplex_ [ i ]; }

  /// operator == 
  bool
  operator == ( Simplex const& rhs ) const;

  friend std::size_t 
  hash_value(Simplex const&);

  friend std::ostream &
  operator << ( std::ostream & outstream,
                Simplex const& simplex );
};

std::ostream & operator << ( std::ostream & outstream, const Simplex & s ) {
  Integer size = s . simplex () . size ();
  outstream << "[";
  for ( Integer i = 0; i < size - 1; ++ i ) {
    outstream << s . simplex () [ i ] << " ";
  }
  if ( s . simplex () . size () > 0 ) outstream << s . simplex () [ size - 1 ]; 
  outstream << "]";
  return outstream;
}

bool Simplex::operator == ( const Simplex & rhs ) const {
  if ( dimension () != rhs . dimension () ) return false;
  for ( unsigned Integer i = 0; i < simplex_ . size (); ++ i ) {
    if ( simplex_ [ i ] != rhs . simplex_ [ i ] ) return false;
  }
  return true;
}

inline
unsigned Integer fvn_hash( unsigned Integer hash_me ) {
  // Fowler Noll Vo hash function.
  unsigned Integer hash = 2166136261u;
  hash ^= (hash_me >> 24);
  hash *= 16777619;
  hash ^= ((hash << 8) >> 24);
  hash *= 16777619;
  hash ^= ((hash_me << 16) >> 24);
  hash *= 16777619;
  hash ^= ((hash_me << 24) >> 24);
  hash *= 16777619;
  return hash;
}

inline
std::size_t hash_value(Simplex const & simplex ) {  
  std::size_t seed = 0;
  Integer d = simplex . dimension ();
  for ( Integer i = 0; i <= d; ++ i ) {
    seed += fvn_hash ( simplex [ i ] );
  }
  return seed;
}


/********************************
 *      SIMPLICIAL COMPLEX      *
 ********************************/
class SimplicialComplex : public Complex {
public:

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

  /// loadFromFile
  ///   the file is a list of simplicies where an n-simplex is a list of Integeregers 
  ///   "Integer_1 Integer_2 ...  Integer_{n+1} \n" each representing a vertex of the complex,
  ///   Assumption: The Integeregers describing each simplex is listed in increasing value
  ///   The function reads the file and adds each simplex as well as all its faces
  ///   No one line can have more than 512 characters
  void
  loadFromFile ( const char * FileName );

  /// loadFromMaxSimplices
  ///   the vector is a list of simplicies where an n-simplex is a list of Integeregers 
  ///   "Integer_1 Integer_2 ...  Integer_{n+1} \n" each representing a vertex of the complex,
  ///   Assumption: The Integeregers describing each simplex is listed in increasing value
  void
  loadFromMaxSimplices (std::vector< std::vector<Integer> > const& max_simplices);
  
private:
  std::unordered_set<Simplex> processed_;
  
  std::vector<Simplex> simplices_;
  std::vector < std::unordered_map<Index, Chain > > coboundaries_;
  
  void generateCoboundaryData ( void );
  
};

/*******************************
 *        DEFINITIONS          *
 *******************************/
inline void SimplicialComplex::
loadFromFile ( const char * FileName) {
  //char *ptr;
  std::ifstream input_file ( FileName ); 
  if ( not input_file . good () ) {
    std::cerr << "SimplicialComplex::loadFromFile. Fatal Error:\n  " << FileName << " not found.\n";
    throw std::runtime_error("File Parsing Error: File not found");
  } /* if */
  //Integer index = 0;
  
  //Simplex s;
  std::vector< std::vector<Integer> > max_simplices;
  std::size_t line_number = 0;
  while ( not input_file . eof () ) {
    ++ line_number;
    std::string line;
    getline( input_file, line );
    // Check that line is sanitized. If not, throw.
    for ( Integer i = 0; i < line.size(); ++ i ) {
      if ( ! ( std::isspace(line[i]) || std::isdigit(line[i]) ) ) {
        std::cerr << "SimplicialComplex::loadFromFile. Fatal Error:\n  Cannot parse line #" << line_number << " of " << FileName << "\n";
        std::cerr << " --> " << line << "\n";
        throw std::runtime_error("File Parsing Error: Invalid file");
      }
    }
    std::vector < Integer > simplex;
    std::istringstream is( line );
    Integer v;
    while ( is >> v ) simplex . push_back ( v );
    std::sort ( simplex . begin (), simplex . end () );
    if ( simplex.size() > 0 ) max_simplices.push_back(simplex);
  }
  input_file . close ();
  
  loadFromMaxSimplices(max_simplices);
}


inline void SimplicialComplex::
loadFromMaxSimplices (const std::vector< std::vector<Integer> > & max_simplices) {
  std::unordered_set < Simplex > processed_;
  std::queue < Simplex > work_q;
  BOOST_FOREACH ( const std::vector < Integer > & maximal, max_simplices ) {
    processed_.insert(maximal);
    if (maximal.size() == 1) continue;
    work_q.push(maximal);
    //std::cout << maximal << "\n";
    while ( not work_q.empty() ) {
      //std::cout << "Front of the queue : " << work_q.front() << "\n";
      for ( unsigned Integer i = 0; i <= work_q.front().dimension(); i ++ ) {
        Simplex work_simplex;
        for (unsigned Integer j = 0; j <= work_q.front().dimension(); j ++ ) {
          if (j != i)
            work_simplex.simplex().push_back(work_q.front()[j]);
        }
        //std::cout << "   Inserting simplex : " << work_simplex << "\n";
        processed_.insert(work_simplex);
        if (work_simplex.dimension() > 0) work_q.push(work_simplex);
      }
      //std::cout << "Popping simplex : " << work_q.front() << "\n";
      work_q.pop();
    }
  }
  
  /*
  BOOST_FOREACH( const Simplex &s, processed_ ) {
    std::cout << s << "\n";
  }
  */
  
  //startInserting ();
  BOOST_FOREACH( const Simplex &s, processed_ ) {
    //std::cout << "Inserting cell " << s << "\n";
    insertCell (s, s.dimension() );
  }
  //finishedInserting ();
  generateCoboundaryData ();
  return;
}


inline void SimplicialComplex::
column ( Integer i, std::function<void(Integer)> const& callback ) const { 
  Simplex s = indexToCell(input, dim);
  //std::cout << "bd(" << input << ", " << dim << ") = ";
  for ( unsigned Integer i = 0; i <= s.dimension(); i++ ) {
    Simplex work_simplex;
    for ( unsigned Integer j = 0; j <= s.dimension(); j++ ) {
      if ( j != i) work_simplex.simplex().push_back( s[j] );
    }
    //std::cout << "Simplex " << work_simplex << "has index " << cellToIndex(work_simplex, dim-1) << "\n";
    (*output) += Term ( cellToIndex (work_simplex, dim-1 ), sign ? positive : negative );
    sign = not sign;
  }
  output -> dimension () = dim - 1;
  //std::cout << *output << "\n";
} /* SimplicialComplex::boundary */

inline void SimplicialComplex::
coboundary ( Chain * output, const Index input, Integer dim ) const {
  //std::cout << "cbd(" << input << ", " << dim << ") = ";

  *output = coboundaries_[dim].find(input)->second;
  //std::cout << *output << "\n";

}

inline void SimplicialComplex::
generateCoboundaryData ( void ) {
// should be called by finalize()
  coboundaries_ . resize ( dimension () + 1 );
  // Insert empty chains.
  for ( Integer dim = 0; dim <= dimension (); ++ dim ) {
    for ( Index i = 0; i < size ( dim ); ++ i ) {
      coboundaries_ [ dim ] . insert ( std::make_pair ( i, Chain () ) );
      coboundaries_ [ dim ] [ i ] . dimension () = dim + 1;
    }
  }
  // Fill in coboundary data
  for ( Integer dim = 0; dim <= dimension (); ++ dim ) {
    for ( Index i = 0; i < size ( dim ); ++ i ) {
      Chain bd = boundary ( i, dim );
      BOOST_FOREACH ( const Term & t, bd () ) {
        coboundaries_ [ dim - 1 ] [ t . index () ] += Term ( i,  t . coef () );
      }
    }
  }
}
  