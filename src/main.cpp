#include <iostream>
#include <Eigen/Dense>

using Eigen::MatrixXd;

// number of basis
int constexpr n = 10;

// distance between basis (nm)
double constexpr base = 0.34;
// distance between elix (nm)
double constexpr sigma_0 = 1.;
// angle between basis (radiants)
double constexpr phi_0 = 0.6283185307179586;
// angle between elix (radiants)
double constexpr omega_0 = 2.217594814298678;


// common length of the rigid rod
double constexpr h_0 = 1.790326582710125;
// equilibrium distance between basis (nm)
double constexpr l_0 = 0.705383591565685;
// equilibrium angle between basis (radiants)
double constexpr theta_0 = 0.5483452644857695;
// coordinates to calculate center coordinates
MatrixXd const d = (MatrixXd(3, 1) << 0, 0, base).finished();

template <typename T>
// coordinates structure
struct Coordinates {
  T x_;
  T y_;
  T z_;

  Coordinates(T x, T y, T z) : x_(x), y_(y), z_(z) {}

  void print() const {
    std::cout << "(" << x_ << " " << y_ << " " << z_ << ")\n";
  }
};

// p1 coordinates
Coordinates<double> const p1 = {sigma_0, 0, 0};
// q1 coordinates
Coordinates<double> const q1 = {sigma_0*std::cos(omega_0), sigma_0*std::sin(omega_0), 0};


template <typename T>
MatrixXd rotation_matrix(T theta, T phi, MatrixXd const& r) {
  T c_phi = std::cos(phi);
  T s_phi = std::sin(phi);
  MatrixXd z_axis(3, 3);
  z_axis << c_phi, -s_phi, 0,
            s_phi, c_phi, 0,
            0, 0, 1;

  T c_theta = std::cos(theta);
  T s_theta = std::sin(theta);
  MatrixXd x_axis(3, 3);
  x_axis << 1, 0, 0,
            0, c_theta, -s_theta,
            0, s_theta, c_theta;
  
  return r * z_axis * x_axis;
}

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
  Base(T theta, T phi, T psi) : theta_(theta), phi_(phi), psi_(psi) {}

  T theta() const { return theta_; }

  T phi() const { return phi_; }

  T psi() const { return psi_; }

  Coordinates<T> central() const { return central_; }

  Coordinates<T> p() const { return p_; }

  Coordinates<T> q() const { return q_; }

  void set_central_coordinates(Coordinates<T> const& coordinates) {
    central_ = coordinates;
  }

  MatrixXd calculate_coordinates(MatrixXd const& previous_rotation, Coordinates<T> const& previous_coordinates) {
    MatrixXd const rotation = rotation_matrix<T>(theta_, phi_, previous_rotation);
    // central coordinates
    MatrixXd const c_m = rotation * d;
    central_.x_ = previous_coordinates.x_ + c_m(0, 0);
    central_.y_ = previous_coordinates.y_ + c_m(1, 0);
    central_.z_ = previous_coordinates.z_ + c_m(2, 0);
    // p coordinates
    MatrixXd const p_m = rotation * (MatrixXd(3, 1) << sigma_0*std::cos(psi_), sigma_0*std::sin(psi_), 0).finished();
    p_.x_ = central_.x_ + p_m(0, 0);
    p_.y_ = central_.y_ + p_m(1, 0);
    p_.z_ = central_.z_ + p_m(2, 0);
    // q coordinates
    MatrixXd const q_m = rotation * (MatrixXd(3, 1) << sigma_0*std::cos(psi_ + omega_0), sigma_0*std::sin(psi_ + omega_0), 0).finished();
    q_.x_ = central_.x_ + q_m(0, 0);
    q_.y_ = central_.y_ + q_m(1, 0);
    q_.z_ = central_.z_ + q_m(2, 0);
    return rotation;
  }
};

int main() {
  Base<double> b_1(0., 0., 0.);
  auto const b_1_c = b_1.central();
  Base<double> b_2(0., 0., 0.);
  {
    // Coordinates<double> const b_2_coo = {0, 0, base};
    // b_2.set_central_coordinates(b_2_coo);
  }
  b_2.calculate_coordinates(MatrixXd::Identity(3, 3), b_1_c);
  auto const b_2_c = b_2.central();
  Base<double> b_3(0., 0., 0.);
  auto const b_3_m = b_3.calculate_coordinates(MatrixXd::Identity(3, 3), b_2_c);
  auto const b_3_c = b_3.central();
  b_1_c.print();
  auto const b_1_p = b_1.p();
  b_1_p.print();
  auto const b_1_q = b_1.q();
  b_1_q.print();
  b_2_c.print();
  auto const b_2_p = b_2.p();
  b_2_p.print();
  auto const b_2_q = b_2.q();
  b_2_q.print();
  b_3_c.print();
  auto const b_3_p = b_3.p();
  b_3_p.print();
  auto const b_3_q = b_3.q();
  b_3_q.print();
}
