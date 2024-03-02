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
double constexpr h_0 = 2 * sigma_0 * std::sin(omega_0 / 2);
// equilibrium distance between basis (nm)
double constexpr l_0 = std::sqrt(2 * sigma_0*sigma_0 * (1 - std::cos(phi_0)) + base*base);
// equilibrium angle between basis (radiants)
double constexpr theta_0 = std::acos((2 * sigma_0*sigma_0 * (1 - std::cos(phi_0)) * std::cos(phi_0) + base*base) / (2 * sigma_0*sigma_0 * (1 - std::cos(phi_0)) + base*base));
// coordinates in eigen
MatrixXd const d = (MatrixXd(3, 1) << 0, 0, base).finished();

template <typename T>
// coordinates structure
struct Coordinates {
  T x_;
  T y_;
  T z_;
};


template <typename T>
MatrixXd rotation_matrix(Angle theta, Angle phi, MatrixXd const& r) {
  T phi = phi.tao_;
  T c_phi = std::cos(phi);
  T s_phi = std::sin(phi);
  MatrixXd z_axis(3, 3);
  z_axis << c_phi, -s_phi, 0,
            s_phi, c_phi, 0,
            0, 0, 1;

  T theta = theta.tao_;
  T c_theta = std::cos(tao);
  T s_theta = std::sin(tao);
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
  Coordinates<T> coordinates_;

  public:
  Base(T theta, T phi, T psi) : theta_(theta), phi_(phi), psi_(psi) {}

  T theta() const { return theta_; }

  T phi() const { return phi_; }

  T psi() const { return psi_; }

  Coordinates<T> coordinates() const { return coordinates_; }

  MatrixXd coordinates(MatrixXd const& r, Coordinates const& c) {
    MatritxXd const rotation = rotation_matrix<T>(theta_, phi_, r);
    coordinates_.x_ = rotation(0, 0) + c.x_;
    coordinates_.y_ = rotation(1, 0) + c.y_;
    coordinates_.z_ = rotation(2, 0) + c.z_;
    return rotation;
  }
};


