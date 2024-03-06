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

  // subtraction of coordinates
  template<typename U>
  Coordinates<T> operator-(Coordinates<U> const& c) const {
    return {x_ - c.x_, y_ - c.y_, z_ - c.z_};
  }

  // scalar division
  template<typename U>
  Coordinates<T> operator/(U const& d) const {
    return {x_ / d, y_ / d, z_ / d};
  }

  // scalar product
  T operator*(Coordinates<T> const& c) const {
    return x_ * c.x_ + y_ * c.y_ + z_ * c.z_;
  }

  // norm of substracted coordinates
  T operator||(Coordinates<T> const& c) const {
    return std::sqrt(std::pow(x_ - c.x_, 2) + std::pow(y_ - c.y_, 2)
                      + std::pow(z_ - c.z_, 2));
  }

  // cross product
  template<typename U>
  Coordinates<T> operator%(Coordinates<U> const& c) const {
    return {y_ * c.z_ - z_ * c.y_, z_ * c.x_ - x_ * c.z_, x_ * c.y_ - y_ * c.x_};
  }

  // addition of coordinates
  template<typename U>
  Coordinates<T>& operator+=(Coordinates<U> const& c) {
    x_ += c.x_;
    y_ += c.y_;
    z_ += c.z_;
    return *this;
  }

  // norm of coordinates
  T norm() const {
    return std::sqrt(x_ * x_ + y_ * y_ + z_ * z_);
  }

  void print() const {
    std::cout << "(" << x_ << " " << y_ << " " << z_ << ")\n";
  }
};
