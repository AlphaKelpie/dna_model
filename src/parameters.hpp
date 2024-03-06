#pragma once

// parameters to create DNA
// vertical distance between basis (nm)
double constexpr base = 0.34;
// distance between elix and backbone structure (nm)
double constexpr sigma_0 = 1.;
// angle between basis (radiants) (36d)
double constexpr psi_0 = 0.6283185307179586;
// angle between P and Q elixes (radiants) (~127.06d)
double constexpr omega_0 = 2.217594814298678;

// parameters to calculate energy
// bonding energy
// constant for bonding energy (pN/nm)
double constexpr k_bond = 1000;
// equilibrium distance between basis of the same elix (nm)
double constexpr l_0 = 0.705383591565685;

// bending energy
// constant for bending energy (pN*nm)
double constexpr k_bend = 7000;
// equilibrium angle between basis of the same elix (radiants) (~31.4d)
double constexpr theta_0 = 0.5483452644857695;

// core energy
// equilibrium distance between the center of the core and each nodal point of the central backbone (nm)
double constexpr sigma_core = 4.5;
// width of the Morse potential well (nm^-1)
double constexpr beta_core = 2.;
// strength of the interaction between the core and the DNA (pNnm)
double constexpr d_core = 10.;

// exclusion energy
// equilibrium distance between not nearby bases (nm)
double constexpr sigma_exc = 2.1;
// width of the Morse potential well exclusion energy (nm^-1)
double constexpr beta_exc = 2.;
// strength of the interaction between not nearby bases (pNnm)
double constexpr d_exc = 1.;
// number of bases to consider as nearby
int constexpr n_nearby = 7;

// physical constants
// Boltzman constant (pN*nm/K)
double constexpr k_b = 1.380649e-2;
// temperature (K)
double constexpr T = 300.;

// parameters to Metropolis algorithm
// step size
double constexpr step = 0.017453292519943295;

// parameter to wraping number
// distance for absorbed bases (nm)
double constexpr absorbed = 4.8;

// parameters to rod energy 1
// strength of the interaction between the rod and the DNA (pNnm)
double constexpr d_rod_1 = 20.;
// width of the Morse potential (nm^-1)
double constexpr beta_rod_1 = 2.;

// parameters to rod energy 2
// strength of the interaction between the rod and the DNA (pNnm)
double constexpr d_rod_2 = 5.;
// width of the Morse potential (nm^-1)
double constexpr beta_rod_2 = 2.;

// parameters to 2 molecules
// parameters to attractive and repulsive energy
// width of the Morse potential well (nm^-1)
double constexpr beta_dna = 0.8;
// equilibrium distance between the two center backbones (nm)
double constexpr sigma_dna = 2.1;
