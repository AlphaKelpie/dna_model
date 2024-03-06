#include "parameters.hpp"
#include "coordinates.hpp"
#include "base.hpp"
#include "functions.hpp"

int main() {
  // number of basis
  int n = 100;
  // Read parameters from file
  std::ifstream file("data.txt");
  if (file.is_open()) {
    file >> n;
    file.close();
  } else {
    std::cerr << "Unable to open file\n";
  }
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
  std::cout << "Energy:\t" << energy << '\n';

  save_coordinates(dna);
}
