#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Dense>

using Eigen::MatrixXd;

// number of basis
int constexpr n = 11;

// vertical distance between basis (nm)
double constexpr base = 0.34;
// distance between elix and backbone structure (nm)
double constexpr sigma_0 = 1.;
// angle between basis (radiants) (36d)
double constexpr psi_0 = 0.6283185307179586;
// angle between P and Q elixes (radiants) (~127.06d)
double constexpr omega_0 = 2.217594814298678;
// constant for bonding energy (pN/nm)
double constexpr k_bond = 1000;
// constant for bending energy (pN*nm)
double constexpr k_bend = 7000;


// common length of the rigid rod between P and Q elixes (nm)
double constexpr h_0 = 1.790326582710125;
// equilibrium distance between basis of the same elix (nm)
double constexpr l_0 = 0.705383591565685;
// equilibrium angle between basis of the same elix (radiants) (~31.4d)
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

// p1 coordinates
Coordinates<double> const p1 = {sigma_0, 0, 0};
// q1 coordinates
Coordinates<double> const q1 = {sigma_0*std::cos(omega_0), sigma_0*std::sin(omega_0), 0};


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

  // void set_central_coordinates(Coordinates<T> const& coordinates) {
  //   central_ = coordinates;
  // }

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

  MatrixXd calculate_coordinates(MatrixXd const& previous_rotation, Coordinates<T> const& previous_coordinates) {
    MatrixXd const rotation = previous_rotation * this->rotation_matrix();
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

template <typename T>
double bonding_energy(Coordinates<T> const& now, Coordinates<T> const& before) {
  return .5 * k_bond * std::pow((now && before) - l_0, 2);
}

template <typename T>
double bending_energy(Coordinates<T> const& now, Coordinates<T> const& before,
                      Coordinates<T> const& after) {
  T const theta = std::acos((now - before) / (now && before)
                            * (after - now) / (after && now));
  return .5 * k_bend * std::pow(theta - theta_0, 2);
}

template <typename T>
double energy(std::vector<Base<T>> const& dna) {
  double energy = 0.;
  for (short int i = 1; i != n; ++i) {
    // bonding energy
    energy = bonding_energy<T>(dna[i].p(), dna[i-1].p())
           + bonding_energy<T>(dna[i].q(), dna[i-1].q());
    // bending energy
    if (i == n-1) { continue; }
    energy += bending_energy<T>(dna[i].p(), dna[i-1].p(), dna[i+1].p())
            + bending_energy<T>(dna[i].q(), dna[i-1].q(), dna[i+1].q());
  }
  return energy;
}

int main() {
  std::vector<Base<double>> dna;
  dna.push_back(Base<double>(0., 0., 0.));
  Coordinates<double> central_previous_coordinates = {0, 0, 0};
  MatrixXd previous_rotation_matrix = MatrixXd::Identity(3, 3);
  for (short int i = 1; i != n; ++i) {
    Base<double> b(0., 0., i*psi_0);
    previous_rotation_matrix = b.calculate_coordinates(previous_rotation_matrix, central_previous_coordinates);
    central_previous_coordinates = b.central();
    dna.push_back(b);
  }
  for (auto const& b : dna) {
    b.central().print();
    b.p().print();
    b.q().print();
    std::cout << "\n\n";
  }
  std::cout << energy(dna) << '\n';
}
