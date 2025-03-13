//
// Created by albpl on 3/12/2025.
//
#include <iostream>

#include "potential_energy_surface.h"
using namespace std;

/// @brief Example of how to use the class PotentialEnergySurface to evaluate
/// the energy of a system
///        with 6 degrees of freedom at a given set of coordinates.
int main() {
  const string path_file_pes1 =
      "../testing_datafiles/coefficients/C1(1)_C1(1)_Coeff.txt";
  PotentialEnergySurface pes1(path_file_pes1);

  const int system_dimension = 6;
  const vector<double> coordinates{12.00, 10, 20, 30, 40, 50};
  const string coordinate_format = "Euler_ZYZ";

  const double energy =
      pes1.EvaluateLRF(system_dimension, coordinates, coordinate_format);
  const double asymptote = pes1.GetAsymptote();
  cout << "Energy: " << energy << endl;
  cout << "Asymptote: " << asymptote << endl;
  return 0;
}
