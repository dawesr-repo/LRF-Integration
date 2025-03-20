//
// Created by albpl on 3/20/2025.
//

#ifndef TENSOR_H
#define TENSOR_H

#include <vector>
#include <iostream>

/// @param dim integer with the number of degrees of freedom of the system
/// @param coordinate_format format of the coordinates(supported formats:
/// "Euler_ZXZ", "Euler_ZYZ" and "Spherical")
/// @param coordinates R, beta1, beta2, alpha, gamma1, gamma2. R must be in
/// Angstroms and the angles in degrees.
class Tensor {
  int system_dimension_ = 6;
  std::vector <double> coordinates_ {10.00,10.0,20.0,30.0,40.0,50.0};
  std::string coordinate_format_ = "Euler_ZYZ";

  public:
  Tensor(const int system_dimension, const std::vector<double> &coordinates,
          const std::string &coordinate_format);
   std::vector<double> UserCoordinatesToGeneralCoordinates() const;
   std::vector<double> CalculateTensor(const int max_t_tensor_order);
   double GetIntermolecularDistance()const {return coordinates_.at(0);}
  };


#endif //TENSOR_H
