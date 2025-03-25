//
// Created by albpl on 3/20/2025.
//

#include <cmath>
#include <iostream>
#include <tuple>
#include <vector>

#include "Tensor.h"

/// @brief ar,br,cab are the projection of the local molecular frame with
/// respect to the laboratory frame after rotation.
/// @return
std::vector<double> Tensor::Ar(const std::vector<double> &general_coordinates) {
  //{Az,Ax,Ay}
  const std::vector<double> &vec_ar = {
      cos(general_coordinates.at(1)),
      sin(general_coordinates.at(1)) * sin(general_coordinates.at(4)),
      cos(general_coordinates.at(4)) * sin(general_coordinates.at(1))};

  return vec_ar;
}

std::vector<double> Tensor::Br(const std::vector<double> &general_coordinates) {
  //{Bz,Bx,By}
  const std::vector<double> &vec_ar = {
      -cos(general_coordinates.at(2)),
      -sin(general_coordinates.at(2)) * sin(general_coordinates.at(5)),
      -cos(general_coordinates.at(5)) * sin(general_coordinates.at(2))};
  return vec_ar;
}

///
///
std::vector<double> Tensor::Cab(
    const std::vector<double> &general_coordinates) {
  //{Az,Ax,Ay}
  const double cos_b1{cos(general_coordinates.at(1))};
  const double sin_b1{sin(general_coordinates.at(1))};
  const double cos_b2{cos(general_coordinates.at(2))};
  const double sin_b2{sin(general_coordinates.at(2))};
  const double cos_phi{cos(general_coordinates.at(3))};
  const double sin_phi{sin(general_coordinates.at(3))};
  const double cos_c1{cos(general_coordinates.at(4))};
  const double sin_c1{sin(general_coordinates.at(4))};
  const double cos_c2{cos(general_coordinates.at(5))};
  const double sin_c2{sin(general_coordinates.at(5))};

  const std::vector<double> &vec_cab = {

      cos_b1 * cos_b2 + cos_phi * sin_b1 * sin_b2,  // Czz
      cos_c2 * sin_phi * sin_b1 +
          (-cos_phi * cos_b2 * sin_b1 + cos_b1 * sin_b2) * sin_c2,  // Czx
      -cos_phi * cos_b2 * cos_c2 * sin_b1 + cos_b1 * cos_c2 * sin_b2 -
          sin_phi * sin_b1 * sin_c2,  // Czy
      cos_b2 * sin_b1 * sin_c1 -
          sin_b2 * (cos_c1 * sin_phi + cos_phi * cos_b1 * sin_c1),  // Cxz
      -cos_b1 * cos_c2 * sin_phi * sin_c1 +
          (cos_b2 * cos_c1 * sin_phi + sin_b1 * sin_b2 * sin_c1) * sin_c2 +
          cos_phi *
              (cos_c1 * cos_c2 + cos_b1 * cos_b2 * sin_c1 * sin_c2),  // Cxx
      cos_c2 * sin_b1 * sin_b2 * sin_c1 +
          cos_b2 * cos_c2 * (cos_c1 * sin_phi + cos_phi * cos_b1 * sin_c1) +
          (-cos_phi * cos_c1 + cos_b1 * sin_phi * sin_c1) * sin_c2,  // Cxy
      cos_b2 * cos_c1 * sin_b1 +
          sin_b2 * (-cos_phi * cos_b1 * cos_c1 + sin_phi * sin_c1),  // Cyz

      cos_c1 * sin_b1 * sin_b2 * sin_c2 +
          cos_b1 * cos_c1 * (-cos_c2 * sin_phi + cos_phi * cos_b2 * sin_c2) -
          sin_c1 * (cos_phi * cos_c2 + cos_b2 * sin_phi * sin_c2),  // Cyx
      -cos_b2 * cos_c2 * sin_phi * sin_c1 +
          cos_c1 * (cos_c2 * sin_b1 * sin_b2 + cos_b1 * sin_phi * sin_c2) +
          cos_phi *
              (cos_b1 * cos_b2 * cos_c1 * cos_c2 + sin_c1 * sin_c2)  // Cyy
  };

  return vec_cab;
}


/// @brief Function that pass from linear index to the real-spherical index
/// @param i :linear index
/// @return component of the T-tensor in real-spherical notation
std::string Tensor::GetSplittingComponent(const int &i) {
  return i < 0 ? "-1" : i == 0 ? "0" : i % 2 == 1 ? "c" : "s";
}

/// @brief Function that pass from real-spherical index to the linear index
int Tensor::GetTensorComponent(const int &mult_ord, const int &k1,
                               const std::string &k2) {
  return (k1 < 0 || mult_ord < 0 || k1 > mult_ord) ? -1
         : k1 == 0                                 ? (k2 == "0" ? 0 : -1)
                   : (k2 == "s" ? 2 * k1 : (k2 == "c" ? 2 * k1 - 1 : 0));
}
/// @brief Auxiliar function for recursive relation function
std::tuple<int, std::string> Tensor::NEta(const std::string &mu, const int &k1,
                                          const std::string &k2) {
  if (k1==0){
    return  std::make_tuple(0, "0");
  }
  if (mu == "x") {
    if (k1 <= 1) {
      return std::make_tuple(0, "0");
    }
    return std::make_tuple(k1-1, k2);

  }
  if (mu == "y"){
    if (k1 <= 1) {
      return std::make_tuple(0, "0");
    }
    if (k2 == "c") {
      return std::make_tuple(k1-1, "s");
    }
    if (k2 == "s") {
      return std::make_tuple(k1-1, "c");
    }
    //return std::make_tuple(k1-1, "0");
  }

  return std::make_tuple(k1, k2);


}
/// @brief Auxiliar function for recursive relation function
double Tensor::Factorial(int n) {
  if (n < 0) {
    return 0.0;
  }
  double res = 1.0;
  for (int i = 2; i <= n; i++) res *= static_cast<double>(i);
  return res;
}
/// @brief Auxiliar function for recursive relation function
double Tensor::FactorialNN(const int &la, const int &ka1, const int &lb,
                           const int &kb1) {
  if (la < 0 || lb < 0 || ka1 < 0 || kb1 < 0 || ka1 > la || kb1 > lb) {
    return 0.0;
  } else {
    return std::sqrt((Factorial(la + ka1) / Factorial(la - ka1)) *
                     (Factorial(lb + kb1) / Factorial(lb - kb1)));
  }
}
/// @brief Auxiliar function for recursive relation function
double Tensor::CoeffM(const std::string &mu, const int &k1,
                      const std::string &k2) {
  double coeff_m = 0.0;

  if (mu == "x") {
    if (k1 == 1) {
      if (k2 == "c") {
        coeff_m = std::sqrt(2.0);
      }
    } else {
      coeff_m = static_cast<double>(k1);
    }
  } else if (mu == "y") {
    if (k1 == 1) {
      if (k2 == "s") {
        coeff_m = std::sqrt(2.0);
      }
    } else {
      if (k2 == "s") {
        coeff_m = static_cast<double>(k1);
      } else {
        coeff_m = -static_cast<double>(k1);
      }
    }
  } else {
    coeff_m = 1.0;
  }

  return coeff_m;
}
