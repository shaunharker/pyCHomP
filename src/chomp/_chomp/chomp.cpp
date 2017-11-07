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
#include "Homology.h"
#include "Fibration.h"
#include "MorseFibration.h"
#include "ConnectionMatrix.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

PYBIND11_MODULE( _chomp, m) {
  ComplexBinding(m);
  CubicalComplexBinding(m);
  MorseMatchingBinding(m);
  MorseComplexBinding(m);
  HomologyBinding(m);
  FibrationBinding(m);
  MorseFibrationBinding(m);
  ConnectionMatrixBinding(m);
}
