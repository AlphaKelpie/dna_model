#pragma once

#include <cmath>
#include <iostream>

template <typename T>
// coordinates structure
struct Coordinates {
  T x_;
  T y_;
  T z_;

  Coordinates(T x, T y, T z) : x_(x), y_(y), z_(z) {}

  T operator-(Coordinates<T> const& c) const {
    return x_ - c.x_ + y_ - c.y_ + z_ - c.z_;
  }

  T operator&&(Coordinates<T> const& c) const {
    return std::sqrt(std::pow(x_ - c.x_, 2) + std::pow(y_ - c.y_, 2)
                      + std::pow(z_ - c.z_, 2));
  }

  void print() const {
    std::cout << "(" << x_ << " " << y_ << " " << z_ << ")\n";
  }
};
