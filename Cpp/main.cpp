//
// Created by albpl on 3/12/2025.
//

#include <iostream>
#include <ostream>
#include <string>

#include "PotentialEnergySurface.h"

int main() {
  const std::string path_file_pes1 =
      "../testing_datafiles/coefficients/C1(1)_C1(1)_Coeff.txt";
  auto* pes1 = new PotentialEnergySurface(path_file_pes1);
  constexpr int system_dimension = 6;
  const std::vector<double> coordinates{10.27, 30.0, 20.0, 120.0, 0.0, 0.0};
  const std::string coordinate_format = "Euler_ZYZ";


  const double energy =
      pes1->EvaluateLRF(system_dimension, coordinates, coordinate_format);

  const double asymptote = pes1->GetAsymptote();

  std::cout << "Energy: " << energy << " (cm^-1)" << std::endl;
  std::cout << "Asymptote: " << asymptote << " (cm^-1)" << std::endl;

  return 0;
}
