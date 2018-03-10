// common.h
/// Shaun Harker 2017-01-25-2341
/// MIT LICENSE

#pragma once

#include <iostream>
#include <cstdint>
#include <numeric>
#include <memory>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <functional>
#include <queue>
#include <limits>
#include <stack>
#include <iterator>
#include "hash.hpp"

// Debug

inline void
print_vector (std::vector<uint64_t> const& v, std::string name) { 
  std::cout << name << " == ["; for ( auto x : v ) std::cout << x << ","; std::cout << "]\n";
}
