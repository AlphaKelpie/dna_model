#pragma once

#include <Eigen/Dense>

#include "parameters.hpp"
#include "coordinates.hpp"

using Eigen::MatrixXd;

// coordinates to calculate center coordinates
MatrixXd const d = (MatrixXd(3, 1) << 0, 0, base).finished();

// p1 coordinates
Coordinates<double> const p1 = {sigma_0, 0, 0};
// q1 coordinates
Coordinates<double> const q1 = {sigma_0*std::cos(omega_0),
                                sigma_0*std::sin(omega_0), 0};

template <typename T>
class Base {
  // bench angle
  T theta_;
  // dihedral angle
  T phi_;
  // twist angle
  T psi_;
  // central coordinates
  Coordinates<T> central_ = {0, 0, 0};
  // p coordinates
  Coordinates<T> p_ = p1;
  // q coordinates
  Coordinates<T> q_ = q1;

  public:
  Base(T theta, T phi, T psi) : theta_{std::fmod(theta, 2*M_PI)},
                                phi_{std::fmod(phi, 2*M_PI)},
                                psi_{std::fmod(psi, 2*M_PI)} {}

  T theta() const { return theta_; }

  T phi() const { return phi_; }

  T psi() const { return psi_; }

  Coordinates<T> central() const { return central_; }

  Coordinates<T> p() const { return p_; }

  Coordinates<T> q() const { return q_; }

  void set_coordinates_c(Coordinates<T> const& c) {
    central_ = c;
  }

  MatrixXd rotation_matrix() {
    T c_phi = std::cos(phi_);
    T s_phi = std::sin(phi_);
    MatrixXd z_axis(3, 3);
    z_axis << c_phi, -s_phi, 0,
              s_phi, c_phi, 0,
              0, 0, 1;

    T c_theta = std::cos(theta_);
    T s_theta = std::sin(theta_);
    MatrixXd x_axis(3, 3);
    x_axis << 1, 0, 0,
              0, c_theta, -s_theta,
              0, s_theta, c_theta;

    return z_axis * x_axis;
  }

  MatrixXd calculate_coordinates(MatrixXd const& previous_rotation,
                                Coordinates<T> const& previous_coordinates) {
    MatrixXd const rotation = previous_rotation * this->rotation_matrix();
    // central coordinates
    MatrixXd const c_m = rotation * d;
    central_.x_ = previous_coordinates.x_ + c_m(0, 0);
    central_.y_ = previous_coordinates.y_ + c_m(1, 0);
    central_.z_ = previous_coordinates.z_ + c_m(2, 0);
    // p coordinates
    MatrixXd const p_m = rotation * (MatrixXd(3, 1) << sigma_0*std::cos(psi_),
                                                       sigma_0*std::sin(psi_),
                                                       0).finished();
    p_.x_ = central_.x_ + p_m(0, 0);
    p_.y_ = central_.y_ + p_m(1, 0);
    p_.z_ = central_.z_ + p_m(2, 0);
    // q coordinates
    MatrixXd const q_m = rotation *
                         (MatrixXd(3, 1) << sigma_0*std::cos(psi_ + omega_0),
                                            sigma_0*std::sin(psi_ + omega_0),
                                            0).finished();
    q_.x_ = central_.x_ + q_m(0, 0);
    q_.y_ = central_.y_ + q_m(1, 0);
    q_.z_ = central_.z_ + q_m(2, 0);
    return rotation;
  }
};
