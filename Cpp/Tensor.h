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
//Helper Functions
  static std::vector<double> Ar(const std::vector<double> &general_coordinates);
  static std::vector<double> Br(const std::vector<double> &general_coordinates);
  static std::vector<double> Cab(const std::vector<double> &general_coordinates);
  static int GetComponent(const int &la, const int &lb, const int &ka, const int &kb);
  static std::string GetSplittingComponent(const int &i);
  static int GetTensorComponent(const int &mult_ord, const int &k1,
                             const std::string &k2);
  static std::tuple<int, std::string> NEta(const std::string &mu, const int &k1, const std::string &k2);
  static double Factorial(int n);
  static double FactorialNN(const int &la, const int &ka1, const int &lb,
                          const int &kb1);
  static double CoeffM(const std::string &mu, const int &k1, const std::string &k2);
  };


#endif //TENSOR_H
