/// MorseMatching.hpp
/// Shaun Harker
/// 2018-02-23
/// MIT LICENSE

#include "MorseMatching.h"
#include "CubicalMorseMatching.h"
#include "GenericMorseMatching.h"

inline
std::shared_ptr<MorseMatching>
MorseMatching::compute_matching ( std::shared_ptr<Complex> complex ) {
  if ( std::dynamic_pointer_cast<CubicalComplex>(complex) ) {
    return std::make_shared<CubicalMorseMatching>(std::dynamic_pointer_cast<CubicalComplex>(complex));
  } else {
    return std::make_shared<GenericMorseMatching>(complex);
  }
}

inline
std::shared_ptr<MorseMatching>
MorseMatching::compute_matching ( std::shared_ptr<Fibration> fibration ) {
  if ( std::dynamic_pointer_cast<CubicalComplex>(fibration->complex()) ) {
    return std::make_shared<CubicalMorseMatching>(fibration);
  } else {
    return std::make_shared<GenericMorseMatching>(fibration);
  }
}
