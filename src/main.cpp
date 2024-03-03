#include <vector>

#include "coordinates.hpp"
#include "base.hpp"

// using Eigen::MatrixXd;

// number of basis
int constexpr n = 11;
// constant for bonding energy (pN/nm)
double constexpr k_bond = 1000;
// constant for bending energy (pN*nm)
double constexpr k_bend = 7000;
// equilibrium distance between basis of the same elix (nm)
double constexpr l_0 = 0.705383591565685;
// equilibrium angle between basis of the same elix (radiants) (~31.4d)
double constexpr theta_0 = 0.5483452644857695;



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
