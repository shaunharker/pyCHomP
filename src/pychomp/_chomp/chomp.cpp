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
#include "GradedComplex.h"
#include "MorseGradedComplex.h"
#include "ConnectionMatrix.h"
#include "Grading.h"
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
  GradedComplexBinding(m);
  MorseGradedComplexBinding(m);
  ConnectionMatrixBinding(m);
  GradingBinding(m);
  SimplicialComplexBinding(m);
  OrderComplexBinding(m);
  DualComplexBinding(m);
}
