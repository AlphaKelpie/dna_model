#pragma once

#include <cmath>
#include <iostream>

template <typename T>
// coordinates structure
struct Coordinates {
  T x_;
  T y_;
  T z_;

  Coordinates() : x_(0), y_(0), z_(0) {}
  Coordinates(T x, T y, T z) : x_(x), y_(y), z_(z) {}

  Coordinates<T> operator-(Coordinates<T> const& c) const {
    return {x_ - c.x_, y_ - c.y_, z_ - c.z_};
  }

  template<typename U>
  Coordinates<T> operator/(U const& d) const {
    return {x_ / d, y_ / d, z_ / d};
  }

  T operator*(Coordinates<T> const& c) const {
    return x_ * c.x_ + y_ * c.y_ + z_ * c.z_;
  }

  T operator&&(Coordinates<T> const& c) const {
    return std::sqrt(std::pow(x_ - c.x_, 2) + std::pow(y_ - c.y_, 2)
                      + std::pow(z_ - c.z_, 2));
  }

  template<typename U>
  Coordinates<T> operator%(Coordinates<U> const& c) const {
    return {y_ * c.z_ - z_ * c.y_, z_ * c.x_ - x_ * c.z_, x_ * c.y_ - y_ * c.x_};
  }

  template<typename U>
  Coordinates<T>& operator+=(Coordinates<U> const& c) {
    x_ += c.x_;
    y_ += c.y_;
    z_ += c.z_;
    return *this;
  }

  void print() const {
    std::cout << "(" << x_ << " " << y_ << " " << z_ << ")\n";
  }
};
