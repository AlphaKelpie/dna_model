#pragma once

#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include <array>
#include <span>

#include "parameters.hpp"
#include "coordinates.hpp"
#include "base.hpp"



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
  int const m = static_cast<int>(dna.size());
  for (int i = 1; i != m; ++i) {
    // bonding energy
    energy += bonding_energy<T>(dna[i].p(), dna[i-1].p())
            + bonding_energy<T>(dna[i].q(), dna[i-1].q());
    // bending energy
    if (i == m-1) { continue; }
    energy += bending_energy<T>(dna[i].p(), dna[i-1].p(), dna[i+1].p())
            + bending_energy<T>(dna[i].q(), dna[i-1].q(), dna[i+1].q());
  }
  return energy;
}

// Save central, p, and q coordinates of Base in separate files
template <typename T>
void save_coordinates(std::vector<Base<T>> const& dna,
                      std::string const& phi = "0",
                      std::string const& theta = "0",
                      std::string const& path = "./") {
  std::ofstream c_file(path + phi + "_" + theta + "_c.txt");
  std::ofstream p_file(path + phi + "_" + theta + "_p.txt");
  std::ofstream q_file(path + phi + "_" + theta + "_q.txt");

  if (c_file.is_open() && p_file.is_open() && q_file.is_open()) {
    c_file << std::fixed << std::setprecision(6);
    p_file << std::fixed << std::setprecision(6);
    q_file << std::fixed << std::setprecision(6);
    for (auto const& b : dna) {
      c_file << b.central().x_ << '\t' << b.central().y_ << '\t'
             << b.central().z_ << '\n';
      p_file << b.p().x_ << '\t' << b.p().y_ << '\t' << b.p().z_ << '\n';
      q_file << b.q().x_ << '\t' << b.q().y_ << '\t' << b.q().z_ << '\n';
    }
    c_file.close();
    p_file.close();
    q_file.close();
  } else {
    std::cerr << "Failed to open coordinates file for writing.\n";
  }
}

// Convert a degree angle to a radian angle
template <typename To, typename From>
To deg2rad(From const& deg) {
  return (To)deg * M_PI / 180;
}

// Save energy in a file
template <typename T, typename U>
void save_energy(std::span<double> const& energies,
                 std::span<T> const& rows = std::array<int, 1>{0},
                 std::span<U> const& columns = std::array<int, 1>{0},
                 std::string const& path = "./") {
  std::ofstream e_file(path + "energy.txt");
  if (e_file.is_open()) {
    e_file << std::fixed << std::setprecision(3);
    int const rows_size = static_cast<int>(rows.size());
    int const columns_size = static_cast<int>(columns.size());
    e_file << "index";
    for (U const& col : columns) {
      e_file << '\t' << col;
    }
    e_file << '\n';
    for (int i = 0; i != rows_size; ++i) {
      e_file << rows[i];
      for (int j = 0; j != columns_size; ++j) {
        e_file << '\t' << energies[i*columns_size + j];
      }
      e_file << '\n';
    }
    e_file.close();
  } else {
    std::cerr << "Failed to open energy file for writing.\n";
  }
}
