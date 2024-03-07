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

int main() {  
  // number of basis
  int n = 100;
  // evolution steps
  int epochs = 2000;
  // path
  std::string path = "./brand/";
  // force on the z axis
  double force = 0.;
  // strength of the interaction
  double d_dna = 0.7;
  // Read parameters from file
  std::ifstream file("data_brand.txt");
  if (file.is_open()) {
    file >> n;
    file >> epochs;
    file >> path;
    file >> force;
    file >> d_dna;
    file.close();
  } else {
    std::cerr << "Default parameters\n";
  }

  path += std::to_string(int(d_dna*100)) + "_"
                          + std::to_string(int(force)) + "_";

  // std::random_device rd;
  // std::mt19937 gen(rd());
  std::mt19937 gen(0);
  std::uniform_int_distribution<int> random_angle(12, 6*(n-1));
  std::uniform_int_distribution<int> oscillation(0, 1);
  std::uniform_real_distribution<double> prob(0., 1.);
  std::cout << std::fixed << std::setprecision(6);

  std::vector<Base<double>> A_dna;
  std::vector<Base<double>> B_dna;
  std::vector<Output> parameters;
  A_dna.push_back(Base<double>(0., 0., 0.));
  B_dna.push_back(Base<double>(0., 0., 0.));
  B_dna[0].set_coordinates(Coordinates<double>(4., 0., 0.));
  A_dna.push_back(Base<double>(0., 0., psi_0));
  B_dna.push_back(Base<double>(0., 0., psi_0));
  Coordinates<double> A_central_previous_coordinates = A_dna[0].central();
  Coordinates<double> B_central_previous_coordinates = B_dna[0].central();
  MatrixXd A_previous_rotation_matrix = MatrixXd::Identity(3, 3);
  MatrixXd B_previous_rotation_matrix = MatrixXd::Identity(3, 3);
  A_dna[1].calculate_coordinates(A_previous_rotation_matrix,
                               A_central_previous_coordinates);
  B_dna[1].calculate_coordinates(B_previous_rotation_matrix,
                                B_central_previous_coordinates);
  A_central_previous_coordinates = A_dna[1].central();
  B_central_previous_coordinates = B_dna[1].central();

  for (int i = 2; i != n; ++i) {
    Base<double> b(0., 0., i*psi_0);
    A_previous_rotation_matrix = b.calculate_coordinates(
                                        A_previous_rotation_matrix,
                                        A_central_previous_coordinates);
    A_central_previous_coordinates = b.central();
    A_dna.push_back(b);
    Base<double> c(0., 0., i*psi_0);
    B_previous_rotation_matrix = c.calculate_coordinates(
                                        B_previous_rotation_matrix,
                                        B_central_previous_coordinates);
    B_central_previous_coordinates = c.central();
    B_dna.push_back(c);
  }
  
  Output params = calculate_branding(A_dna, B_dna, force, d_dna);
  parameters.push_back(params);

  for (int e : tq::trange(epochs)) {
    if (e % 1000 == 0) {
      save_coordinates(A_dna, path + "A_" + std::to_string(e) + "_");
      save_coordinates(B_dna, path + "B_" + std::to_string(e) + "_");
    }
    std::vector<Base<double>> A_dna_new{A_dna[0], A_dna[1]};
    std::vector<Base<double>> B_dna_new{B_dna[0], B_dna[1]};
    A_previous_rotation_matrix = MatrixXd::Identity(3, 3);
    B_previous_rotation_matrix = MatrixXd::Identity(3, 3);
    A_central_previous_coordinates = A_dna_new[1].central();
    B_central_previous_coordinates = B_dna_new[1].central();
    int const random_num = random_angle(gen);
    int const index_base_change = random_num / 6;
    for (int i = 2; i != n; ++i) {
      if (i != index_base_change) {
        Base<double> b(A_dna[i].theta(), A_dna[i].phi(), A_dna[i].psi());
        A_previous_rotation_matrix = b.calculate_coordinates(
                                            A_previous_rotation_matrix,
                                            A_central_previous_coordinates);
        A_central_previous_coordinates = b.central();
        A_dna_new.push_back(b);
        Base<double> c(B_dna[i].theta(), B_dna[i].phi(), B_dna[i].psi());
        B_previous_rotation_matrix = c.calculate_coordinates(
                                            B_previous_rotation_matrix,
                                            B_central_previous_coordinates);
        B_central_previous_coordinates = c.central();
        B_dna_new.push_back(c);
      } else {
        auto A_theta = A_dna[i].theta();
        auto A_phi = A_dna[i].theta();
        auto A_psi = A_dna[i].theta();
        auto B_theta = B_dna[i].theta();
        auto B_phi = B_dna[i].theta();
        auto B_psi = B_dna[i].theta();
        switch (random_num % 6)
        {
          case 0:
            A_theta += (oscillation(gen)*2 - 1) * static_cast<double>(step);
            break;
          case 1:
            A_phi += (oscillation(gen)*2 - 1) * static_cast<double>(step);
            break;
          case 2:
            A_psi += (oscillation(gen)*2 - 1) * static_cast<double>(step);
            break;
          case 3:
            B_theta += (oscillation(gen)*2 - 1) * static_cast<double>(step);
            break;
          case 4:
            B_phi += (oscillation(gen)*2 - 1) * static_cast<double>(step);
            break;
          case 5:
            B_psi += (oscillation(gen)*2 - 1) * static_cast<double>(step);
            break;
          default:
            break;  // gestire eccezione
        }
        Base<double> b(A_theta, A_phi, A_psi);
        A_previous_rotation_matrix = b.calculate_coordinates(
                                            A_previous_rotation_matrix,
                                            A_central_previous_coordinates);
        A_central_previous_coordinates = b.central();
        A_dna_new.push_back(b);
        Base<double> c(B_theta, B_phi, B_psi);
        B_previous_rotation_matrix = c.calculate_coordinates(
                                            B_previous_rotation_matrix,
                                            B_central_previous_coordinates);
        B_central_previous_coordinates = c.central();
        B_dna_new.push_back(c);
      }
    }
    Output const params_new = calculate_branding(A_dna_new, B_dna_new, force,
                                                 d_dna);
    parameters.push_back(params_new);
    if (p(params_new.e_, params.e_) > prob(gen)) {
      A_dna = std::move(A_dna_new);
      B_dna = std::move(B_dna_new);
      params = params_new;
    }
  }

  save_coordinates(A_dna, path + "A_" + std::to_string(epochs) + "_");
  save_coordinates(B_dna, path + "B_" + std::to_string(epochs) + "_");
  std::vector<int> cicles(epochs);
  std::iota(cicles.begin(), cicles.end(), 0);
  std::array<int, 1> useless = {0};
  save_parameters<int, int>(parameters, cicles, useless, path);
  std::cout << '\n';
}
