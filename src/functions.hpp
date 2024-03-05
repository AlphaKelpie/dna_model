#pragma once

#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include <array>
#include <span>

#include "parameters.hpp"
#include "coordinates.hpp"
#include "output.hpp"
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
double core_energy(Coordinates<T> const& base, Coordinates<T> const& core) {
  double const dist = double(base && core) - sigma_core;
  return d_core * (std::exp(-2. * beta_core * dist) - 2. * std::exp(-beta_core * dist));
}

template <typename T>
double exclusion_energy(Coordinates<T> const& base, std::vector<Base<T>> const& over) {
  double energy = 0.;
  for (auto const& b : over) {
    double const dist = double(base && b.central()) - sigma_exc;
    energy += d_exc * std::exp(-2. * beta_exc * dist);
  }
  return energy;
}

template <typename T>
int nearby(Coordinates<T> const& base, Coordinates<T> const& core) {
  return (base && core) <= absorbed;
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

template <typename T>
Output calculate_parameters(std::vector<Base<T>> const& dna, Coordinates<T> const& core) {
  double energy = 0.;
  int near = 0;
  int const m = static_cast<int>(dna.size());
  // core energy for the first base
  energy += core_energy<T>(dna[0].central(), core);
  std::vector<Base<T>> over(dna.begin() + 1 + n_nearby, dna.end());
  // exclusion energy for the first base
  energy += exclusion_energy<T>(dna[0].central(), over);
  // if the first base is nearby
  near += nearby<T>(dna[0].central(), core);
  for (int i{1}; i != m; ++i) {
    // bonding energy
    energy += bonding_energy<T>(dna[i].p(), dna[(i-1)].p())
            + bonding_energy<T>(dna[i].q(), dna[(i-1)].q());
    // core energy
    energy += core_energy<T>(dna[i].central(), core);
    // if the base is nearby
    near += nearby<T>(dna[i].central(), core);
    if (i >= m-1) { continue; }
    // bending energy
    energy += bending_energy<T>(dna[i].p(), dna[i-1].p(), dna[i+1].p())
            + bending_energy<T>(dna[i].q(), dna[i-1].q(), dna[i+1].q());
    if (i >= m-n_nearby-1) { continue; }
    // exclusion energy
    over.erase(over.begin());
    energy += exclusion_energy<T>(dna[i].central(), over);
  }
  double const wrapping_num = b * (near - 1) / (2. * M_PI * sigma_core);
  return {energy, wrapping_num, 0.};
}

// Save central, p, and q coordinates of Base in separate files
template <typename T>
void save_coordinates(std::vector<Base<T>> const& dna,
                      std::string const& path = "./",
                      std::string const& phi = "0",
                      std::string const& theta = "0") {
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
void save_parameters(std::span<Output> const& params,
                 std::span<T> const rows,
                 std::span<U> const columns,
                 std::string const& path = "./") {
  std::ofstream e_file(path + "energy.txt");
  std::ofstream w_file(path + "wrapping.txt");
  std::ofstream c_file(path + "chirality.txt");
  if (e_file.is_open() && w_file.is_open() && c_file.is_open()) {
    e_file << std::fixed << std::setprecision(3);
    w_file << std::fixed << std::setprecision(3);
    c_file << std::fixed << std::setprecision(3);
    int const rows_size = static_cast<int>(rows.size());
    int const columns_size = static_cast<int>(columns.size());
    e_file << "index";
    w_file << "index";
    c_file << "index";
    for (U const& col : columns) {
      e_file << '\t' << col;
      w_file << '\t' << col;
      c_file << '\t' << col;
    }
    e_file << '\n';
    w_file << '\n';
    c_file << '\n';
    for (int i = 0; i != rows_size; ++i) {
      e_file << rows[i];
      w_file << rows[i];
      c_file << rows[i];
      for (int j = 0; j != columns_size; ++j) {
        e_file << '\t' << params[i*columns_size + j].e_;
        w_file << '\t' << params[i*columns_size + j].w_;
        c_file << '\t' << params[i*columns_size + j].c_;
      }
      e_file << '\n';
      w_file << '\n';
      c_file << '\n';
    }
    e_file.close();
    w_file.close();
    c_file.close();
  } else {
    std::cerr << "Failed to open parameters files for writing.\n";
  }
};
