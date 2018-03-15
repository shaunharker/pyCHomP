/// chompy.cpp
/// Shaun Harker
/// 2017-07-20
/// MIT LICENSE

#include "Integer.h"
#include "Iterator.h"
#include "Chain.h"
#include "Complex.h"
#include "CubicalComplex.h"
#include "MorseComplex.h"
#include "MorseMatching.h"
#include "MorseMatching.hpp"
#include "CubicalMorseMatching.h"
#include "GenericMorseMatching.h"
#include "Homology.h"
#include "Fibration.h"
#include "MorseFibration.h"
#include "ConnectionMatrix.h"
#include "Valuation.h"
#include "SimplicialComplex.h"
#include "OrderComplex.h"
#include "DualComplex.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

PYBIND11_MODULE( _chomp, m) {
  ComplexBinding(m);
  CubicalComplexBinding(m);
  MorseMatchingBinding(m);
  CubicalMorseMatchingBinding(m);
  GenericMorseMatchingBinding(m);
  MorseComplexBinding(m);
  HomologyBinding(m);
  FibrationBinding(m);
  MorseFibrationBinding(m);
  ConnectionMatrixBinding(m);
  ValuationBinding(m);
  SimplicialComplexBinding(m);
  OrderComplexBinding(m);
  DualComplexBinding(m);
}
