#include <vector>
#include <fstream>
#include <iomanip>
#include <random>
#include <numeric>

#include "tqdm.hpp"
#include "coordinates.hpp"
#include "output.hpp"
#include "base.hpp"
#include "functions.hpp"

// evolution steps
int constexpr epochs = 100000;

// number of basis
int constexpr n = 200;

// coordinates of the histone (nm)
Coordinates<double> const histone = {0., -4.5, 10.2};

// path
std::string const path = "./histone/";

int main() {
  // std::random_device rd;
  // std::mt19937 gen(rd());
  std::mt19937 gen(0);
  std::uniform_int_distribution<int> random_angle(6, 3*(n-1));
  std::uniform_int_distribution<int> oscillation(0, 1);
  std::uniform_real_distribution<double> prob(0., 1.);
  std::cout << std::fixed << std::setprecision(6);

  std::vector<Base<double>> dna;
  std::vector<Output> parameters;
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
  
  Output params = calculate_parameters(dna, histone);
  parameters.push_back(params);

  for (int e : tq::trange(epochs)) {
    if (e % 1000 == 0) {
      save_coordinates(dna, path + std::to_string(e) + "_");
    }
    std::vector<Base<double>> dna_new{dna[0], dna[1]};
    previous_rotation_matrix = MatrixXd::Identity(3, 3);
    central_previous_coordinates = dna_new[1].central();
    int const random_num = random_angle(gen);
    int const index_base_change = random_num / 3;
    for (int i = 2; i != n; ++i) {
      if (i != index_base_change) {
        Base<double> b(dna[i].theta(), dna[i].phi(), dna[i].psi());
        previous_rotation_matrix = b.calculate_coordinates(previous_rotation_matrix,
                                        central_previous_coordinates);
        central_previous_coordinates = b.central();
        dna_new.push_back(b);
      } else {
        auto theta = dna[i].theta();
        auto phi = dna[i].phi();
        auto psi = dna[i].psi();
        switch (random_num % 3)
        {
          case 0:
            theta += (oscillation(gen)*2 - 1) * static_cast<double>(step);
            break;
          case 1:
            phi += (oscillation(gen)*2 - 1) * static_cast<double>(step);
            break;
          case 2:
            psi += (oscillation(gen)*2 - 1) * static_cast<double>(step);
            break;
          default:
            break;  // gestire eccezione
        }
        Base<double> b(theta, phi, psi);
        previous_rotation_matrix = b.calculate_coordinates(previous_rotation_matrix,
                                        central_previous_coordinates);
        central_previous_coordinates = b.central();
        dna_new.push_back(b);
      }
    }
    Output const params_new = calculate_parameters(dna_new, histone);
    parameters.push_back(params_new);
    if (p(params_new.e_, params.e_) > prob(gen)) {
      dna = std::move(dna_new);
      params = params_new;
    }
  }

  save_coordinates(dna, path + std::to_string(epochs) + "_");
  std::vector<int> cicles(epochs);
  std::iota(cicles.begin(), cicles.end(), 0);
  std::array<int, 1> useless = {0};
  save_parameters<int, int>(parameters, cicles, useless, path);
  std::cout << '\n';
}
