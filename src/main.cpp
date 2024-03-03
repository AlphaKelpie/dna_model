#include <vector>
#include <fstream>
#include <iomanip>
#include <random>

#include "tqdm.hpp"
#include "coordinates.hpp"
#include "base.hpp"

// using Eigen::MatrixXd;

// evolution steps
int constexpr epochs = 100000;
// step size
double constexpr step = 0.017453292519943295;

// number of basis
int constexpr n = 100;
// constant for bonding energy (pN/nm)
double constexpr k_bond = 1000;
// constant for bending energy (pN*nm)
double constexpr k_bend = 7000;
// equilibrium distance between basis of the same elix (nm)
double constexpr l_0 = 0.705383591565685;
// equilibrium angle between basis of the same elix (radiants) (~31.4d)
double constexpr theta_0 = 0.5483452644857695;
// Boltzman constant (pN*nm/K)
double constexpr k_b = 1.380649e-2;
// temperature (K)
double constexpr T = 300.;

// probability function
double p(double const& now, double const& before) {
  return std::exp(-(now - before) / (k_b * T));
}

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
double calculate_energy(std::vector<Base<T>> const& dna) {
  double energy = 0.;
  for (int i = 1; i != n; ++i) {
    // bonding energy
    energy += bonding_energy<T>(dna[i].p(), dna[i-1].p())
            + bonding_energy<T>(dna[i].q(), dna[i-1].q());
    // bending energy
    if (i == n-1) { continue; }
    energy += bending_energy<T>(dna[i].p(), dna[i-1].p(), dna[i+1].p())
            + bending_energy<T>(dna[i].q(), dna[i-1].q(), dna[i+1].q());
  }
  return energy;
}

int main() {
  // std::random_device rd;
  // std::mt19937 gen(rd());
  std::mt19937 gen(0);
  std::uniform_int_distribution<int> oscillation(-1, 1);
  std::uniform_real_distribution<double> uniform(0., 1.);
  std::cout << std::fixed << std::setprecision(6);
  std::vector<Base<double>> dna;
  dna.push_back(Base<double>(0., 0., 0.));
  dna.push_back(Base<double>(0., 0., psi_0));
  Coordinates<double> central_previous_coordinates = dna[0].central();
  MatrixXd previous_rotation_matrix = MatrixXd::Identity(3, 3);
  dna[1].calculate_coordinates(previous_rotation_matrix,
                               central_previous_coordinates);
  central_previous_coordinates = dna[1].central();

  for (int i = 2; i != n; ++i) {
    Base<double> b(0., 0., i*psi_0);
    previous_rotation_matrix = b.calculate_coordinates(previous_rotation_matrix,
                                    central_previous_coordinates);
    central_previous_coordinates = b.central();
    dna.push_back(b);
  }
  double energy = calculate_energy(dna);

  for (int e : tq::trange(epochs)) {
    std::vector<Base<double>> dna_new{dna[0], dna[1]};
    previous_rotation_matrix = MatrixXd::Identity(3, 3);
    central_previous_coordinates = dna_new[1].central();
    for (int i = 2; i != n; ++i) {
      double const theta = dna[i].theta() + oscillation(gen) * step;
      double const phi = dna[i].phi() + oscillation(gen) * step;
      double const psi = dna[i].psi() + oscillation(gen) * step;
      Base<double> b(theta, phi, psi);
      previous_rotation_matrix = b.calculate_coordinates(previous_rotation_matrix,
                                      central_previous_coordinates);
      central_previous_coordinates = b.central();
      dna_new.push_back(b);
    }
    double const energy_new = calculate_energy(dna_new);
    if (p(energy_new, energy) > uniform(gen)) {
      dna = std::move(dna_new);
      energy = energy_new;
    }
  }

  // Save central, p, and q coordinates of Base in separate files
  std::ofstream c_file("./file_c.txt");
  std::ofstream p_file("./file_p.txt");
  std::ofstream q_file("./file_q.txt");

  if (c_file.is_open() && p_file.is_open() && q_file.is_open()) {
    c_file << std::fixed << std::setprecision(6);
    p_file << std::fixed << std::setprecision(6);
    q_file << std::fixed << std::setprecision(6);
    for (auto const& b : dna) {
      c_file << b.central().x_ << '\t' << b.central().y_ << '\t' << b.central().z_ << '\n';
      p_file << b.p().x_ << '\t' << b.p().y_ << '\t' << b.p().z_ << '\n';
      q_file << b.q().x_ << '\t' << b.q().y_ << '\t' << b.q().z_ << '\n';
    }
    c_file.close();
    p_file.close();
    q_file.close();
  } else {
    std::cerr << "Failed to open file for writing.\n";
  }
}
