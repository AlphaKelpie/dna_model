#include <algorithm>

#include "tqdm.hpp"
#include "parameters.hpp"
#include "coordinates.hpp"
#include "base.hpp"
#include "functions.hpp"

int main(int argc, char* argv[]) {
  // number of basis
  int n = 100;
  // path
  std::string path = "./bend_writhe/";

  // Read parameters from file
  if (argc == 2) {
    std::ifstream file("data_bend.txt");
    if (file.is_open()) {
      file >> n;
      file >> path;
      file.close();
    } else {
      std::cerr << "Default parameters\n";
    }
  } else {
    std::cerr << "Default parameters\n";
  }
  
  std::cout << std::fixed << std::setprecision(6);
  
  // simulations for dna graphics
  std::array<int, 2> phies_a = {-4, 4};
  std::array<int, 4> thetas = {0, 1, 3, 5};
  for (int phi : phies_a) {
    for (int theta : thetas) {
      std::vector<Base<double>> dna;
      dna.push_back(Base<double>(0., 0., 0.));
      dna.push_back(Base<double>(0., 0., static_cast<double>(psi_0)));
      auto central_previous_coordinates = dna[0].central();
      MatrixXd previous_rotation_matrix = MatrixXd::Identity(3, 3);
      dna[1].calculate_coordinates(previous_rotation_matrix,
                                   central_previous_coordinates);
      central_previous_coordinates = dna[1].central();
      for (int i = 2; i != n; ++i) {
        Base<double> b(deg2rad<double>(theta), deg2rad<double>(phi),
                       i*static_cast<double>(psi_0));
        previous_rotation_matrix = b.calculate_coordinates(
                                        previous_rotation_matrix,
                                        central_previous_coordinates);
        central_previous_coordinates = b.central();
        dna.push_back(b);
      }

      save_coordinates(dna, path, std::to_string(phi),
                       std::to_string(theta));
    }
  }

  // simulations for energy graphics
  std::vector<double> energies;
  std::vector<double> phies_v(8001);
  double seed = -4.000;
  double const step = 0.001;
  std::generate(phies_v.begin(), phies_v.end(), [&seed, step]() {
    double ret = seed;
    seed += step;
    return ret;
  });

  for (double phi : tq::tqdm(phies_v)) {
    for (int theta : thetas) {
      std::vector<Base<double>> dna;
      dna.push_back(Base<double>(0., 0., 0.));
      dna.push_back(Base<double>(0., 0., static_cast<double>(psi_0)));
      auto central_previous_coordinates = dna[0].central();
      MatrixXd previous_rotation_matrix = MatrixXd::Identity(3, 3);
      dna[1].calculate_coordinates(previous_rotation_matrix,
                                   central_previous_coordinates);
      central_previous_coordinates = dna[1].central();
      for (int i = 2; i != n; ++i) {
        Base<double> b(deg2rad<double>(theta), deg2rad<double>(phi),
                       i*static_cast<double>(psi_0));
        previous_rotation_matrix = b.calculate_coordinates(
                                        previous_rotation_matrix,
                                        central_previous_coordinates);
        central_previous_coordinates = b.central();
        dna.push_back(b);
      }
      energies.push_back(calculate_energy(dna));
    }
  }

  save_energy<double, int>(energies, phies_v, thetas, path);
  std::cout << '\n';
}
