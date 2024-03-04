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
// constant for bonding energy (pN/nm)
double constexpr k_bond = 1000;
// equilibrium distance between basis of the same elix (nm)
double constexpr l_0 = 0.705383591565685;
// constant for bending energy (pN*nm)
double constexpr k_bend = 7000;
// equilibrium angle between basis of the same elix (radiants) (~31.4d)
double constexpr theta_0 = 0.5483452644857695;

// physical constants
// Boltzman constant (pN*nm/K)
double constexpr k_b = 1.380649e-2;
// temperature (K)
double constexpr T = 300.;

// parameters to Metropolis algorithm
// step size
double constexpr step = 0.017453292519943295;
