/// Iterator.h
/// Shaun Harker
/// 2017-07-19
/// MIT LICENSE

#pragma once

#include <boost/iterator/counting_iterator.hpp>
#include <boost/range.hpp>

typedef boost::counting_iterator<Integer> Iterator;
typedef boost::iterator_range<Iterator> Range;
