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
// coordinates in eigen
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
  
  return r * z_axis * x_axis * d;
}

template <typename T>
class Base {
  // bench angle
  T theta_;
  // dihedral angle
  T phi_;
  // twist angle
  T psi_;
  // coordinates
  Coordinates<T> coordinates_ = {0, 0, 0};

  public:
  Base(T theta, T phi, T psi) : theta_(theta), phi_(phi), psi_(psi) {}

  T theta() const { return theta_; }

  T phi() const { return phi_; }

  T psi() const { return psi_; }

  Coordinates<T> coordinates() const { return coordinates_; }

  void set_coordinates(Coordinates<T> const& c) {
    coordinates_ = c;
  }

  MatrixXd calculate_coordinates(MatrixXd const& r, Coordinates<T> const& c) {
    MatrixXd const rotation = rotation_matrix<T>(theta_, phi_, r);
    coordinates_.x_ = rotation(0, 0) + c.x_;
    coordinates_.y_ = rotation(1, 0) + c.y_;
    coordinates_.z_ = rotation(2, 0) + c.z_;
    return rotation;
  }
};

int main() {
  Base<double> b_1(0., 0., 0.);
  auto const b_1_coo = b_1.coordinates();
  Base<double> b_2(0., 0., 0.);
  {
    Coordinates<double> const b_2_coo = {0, 0, base};
    b_2.set_coordinates(b_2_coo);
  }
  auto const b_2_coo = b_2.coordinates();
  Base<double> b_3(0.3, 0., 0.);
  auto const b_3_m = b_3.calculate_coordinates(MatrixXd::Identity(3, 3), b_2_coo);
  auto const b_3_coo = b_3.coordinates();
  b_1_coo.print();
  b_2_coo.print();
  b_3_coo.print();
}
