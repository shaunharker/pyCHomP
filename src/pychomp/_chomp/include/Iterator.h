/// Iterator.h
/// Shaun Harker
/// 2018-03-10
/// MIT LICENSE

#pragma once

#include "common.h"

#include "Integer.h"

template < typename I >
class IteratorRange {
public:
  typedef I iterator;
  typedef I const_iterator;
  IteratorRange ( void ) {}
  IteratorRange ( iterator b, iterator e ) : begin_(b), end_(e) {}
  iterator begin ( void ) const { return begin_;}
  iterator end ( void ) const { return end_;}
  uint64_t size ( void ) const { return end_ - begin_;}
  typename iterator::value_type operator [] ( int64_t i ) const {return *(begin_ + i);}
private:
  iterator begin_;
  iterator end_;
};

class CountingIterator {
public:
    typedef CountingIterator self_type;
    typedef Integer value_type;
    typedef Integer& reference;
    typedef Integer* pointer;
    typedef Integer difference_type;
    typedef std::forward_iterator_tag iterator_category;
    CountingIterator(void) : val_(0) {}
    CountingIterator(Integer i) : val_(i) { }
    self_type operator++() { return ++ val_; }
    self_type operator++(int) { return val_ ++;  }
    self_type operator+(Integer i) const {return CountingIterator(val_ + i);}
    difference_type operator-(self_type const& rhs) const{return val_ - rhs.val_;}
    Integer operator*() const { return val_; }
    //const T* operator->() { return ptr_; }
    self_type operator=(const self_type& other) { val_ = other.val_; return *this; }
    bool operator==(const self_type& rhs)const { return val_ == rhs.val_; }
    bool operator!=(const self_type& rhs)const { return val_ != rhs.val_; }
private:
    int64_t val_;
};

typedef CountingIterator Iterator;
typedef IteratorRange<Iterator> Range;
