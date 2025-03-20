//
// Created by albpl on 3/20/2025.
//

#include "Tensor.h"

#include <cfloat>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <tuple>
#include <vector>

using Vector = std::vector<double>;
//******************* Tensor Methods *******************************************

/// @param system_dimension integer with the number of degrees of freedom of the
/// system
/// @param coordinate_format format of the coordinates(supported formats:
/// "Euler_ZXZ", "Euler_ZYZ" and "Spherical")
/// @param coordinates R, beta1, beta2, alpha, gamma1, gamma2. R must be in
/// Angstroms and the angles in degrees.
Tensor::Tensor(const int system_dimension, const Vector &coordinates,
               const std::string &coordinate_format) {
  system_dimension_ = system_dimension;
  coordinates_ = coordinates;
  coordinate_format_ = coordinate_format;
}

/// @brief Pass from user coordinates to general coordinates
/// @return the standard 6D vector in Euler_ZXZ vector.
///
Vector Tensor::UserCoordinatesToGeneralCoordinates() const {
  // general_coordinates
  Vector vec(6, 0);
  if (system_dimension_ >= 2 && system_dimension_ <= 6) {
    vec.at(0) = coordinates_.at(0);  // R
    vec.at(1) = coordinates_.at(1);  // beta1
    if (system_dimension_ == 3) {
      vec.at(4) = coordinates_.at(2);  // gamma1
    } else if (system_dimension_ >= 4) {
      vec.at(2) = coordinates_.at(2);  // beta2
      vec.at(3) = coordinates_.at(3);  // alpha
    }
    if (system_dimension_ == 5) {
      vec.at(4) = coordinates_.at(4);  // gamma1
    }
    if (system_dimension_ == 6) {
      vec.at(5) = coordinates_.at(5);  // gamma2
    }
  } else {
    std::cout << "Wrong dimension: " << system_dimension_ << std::endl;
    // throw 1;
  }

  if (coordinate_format_ == "Euler_ZYZ") {
    vec.at(4) = vec.at(4) - 90.0;
    vec.at(5) = vec.at(5) - 90.0;
  } else if (coordinate_format_ == "Spherical") {
    vec.at(4) = 90.0 - vec.at(4);
    vec.at(5) = 90.0 - vec.at(5);
  }
  const double pii = acos(-1.0);
  for (int i = 1; i < 6; ++i) {
    vec.at(i) = vec.at(i) * (pii / 180.0);
  }

  return vec;
}

/// @brief Calculate the T-tensor components for a given set of coordinates
/// @param max_t_tensor_order: Maximum order to calculate the T-Tensors. It must
///                            be a positive integer
///
/// @return 1D-vector with the T-tensor components (4-D tensor reshaped)

Vector Tensor::CalculateTensor(const int max_t_tensor_order) {
  std::vector<std::string> coord{"z", "x", "y"};  //! Cartesian Axis Labels

  if (coordinates_.size() != system_dimension_) {
    std::cout << "coordinates size must be equal to dimension" << std::endl;
    // throw 0;
  }
  if (coordinates_.at(0) < 1) {
    std::cout << "coordinates size must be equal to dimension" << std::endl;
    // throw 0;
  }

  const auto general_coordinates_zxz = UserCoordinatesToGeneralCoordinates();

  auto a = Ar(general_coordinates_zxz);
  auto b = Br(general_coordinates_zxz);
  auto cc = Cab(general_coordinates_zxz);

  const int cpns_per_order[15]{1,    7,    26,   70,   155,  301,  532, 876,
                               1365, 2035, 2926, 4082, 5551, 7385, 9640};

  Vector t_tensor(cpns_per_order[max_t_tensor_order - 1], 0);

  int cpn{0};

  for (int order = 1; order <= max_t_tensor_order; ++order) {
    for (int la = 0; la <= order - 1; ++la) {
      int lb = order - la - 1;
      for (int ka = 0; ka <= 2 * la; ++ka) {
        for (int kb = 0; kb <= 2 * lb; ++kb) {
          // Calculating T-Tensor Component
          const std::string ka2 = GetSplittingComponent(ka);
          const std::string kb2 = GetSplittingComponent(kb);

          const int ka1 = floor((ka + 1.0) / 2.0);
          const int kb1 = floor((kb + 1.0) / 2.0);

          if (la == 0 && lb == 0) {  // 1st order tensor component
            t_tensor.at(0) = 1.0;
          } else {
            // recursive relation for lb = 0
            if (lb == 0) {
              // initializing component
              double comp_lk{0};
              const auto la_fact = static_cast<double>((2.0 * la - 1.0) / la);

              // loop though every coordinate axis
              for (int i = 1; i <= 3; i++) {
                // new multipole components
                int rk1;
                std::string rk2;
                tie(rk1, rk2) = NEta(coord.at(i - 1), ka1, ka2);
                const int rk_ = GetTensorComponent(la, rk1, rk2);
                const double m = CoeffM(coord.at(i - 1), ka1, ka2);
                // coefficient NN of the recurrence
                const double fact_nn = FactorialNN(la - 1, rk1, 0, 0);

                if (std::abs(m) > DBL_EPSILON && la >= 1 &&
                    fact_nn > DBL_EPSILON && rk_ <= 2 * (la - 1)) {
                  const int t_cpn = GetComponent(la - 1, 0, rk_, 0);
                  double prod_comp = a.at(i - 1) * t_tensor.at(t_cpn);
                  double fact_prod = la_fact * m * fact_nn;
                  comp_lk = comp_lk + fact_prod * prod_comp;
                }
              }
              // Second term of the recurrence
              if (la >= 2 && ka <= 2 * (la - 2) && ka >= 0) {
                const auto la2_fact = static_cast<double>((la - 1.0) / la);

                const int t_cpn = GetComponent(la - 2, 0, ka, 0);
                comp_lk = comp_lk - la2_fact * FactorialNN(la - 2, ka1, 0, 0) *
                                        t_tensor.at(t_cpn);
              }

              t_tensor.at(cpn) = comp_lk / FactorialNN(la, ka1, lb, kb1);
            }
            // recursive relation for la = 0
            else if (la == 0) {
              // initializing component
              double comp_lk{0};
              const auto lb_fact = static_cast<double>((2.0 * lb - 1.0) / lb);

              // loop though every coordinate axis
              for (int i = 1; i <= 3; i++) {
                // new multipole components
                int rk1;
                std::string rk2;
                tie(rk1, rk2) = NEta(coord.at(i - 1), kb1, kb2);
                const int rk_ = GetTensorComponent(lb, rk1, rk2);
                const double m = CoeffM(coord.at(i - 1), kb1, kb2);
                // coefficient NN of the recurrence
                const double fact_nn = FactorialNN(0, 0, lb - 1, rk1);

                if (std::abs(m) > DBL_EPSILON && lb >= 1 &&
                    fact_nn > DBL_EPSILON && rk_ <= 2 * (lb - 1) && rk_ >= 0) {
                  const int t_cpn = GetComponent(0, lb - 1, 0, rk_);
                  comp_lk = comp_lk + lb_fact * m * fact_nn * b.at(i - 1) *
                                          t_tensor.at(t_cpn);
                }
              }
              // Second term of the recurrence
              if (lb >= 2 && kb <= 2 * (lb - 2) && kb >= 0) {
                const auto lb2_fact = static_cast<double>((lb - 1.0) / lb);
                const int t_cpn = GetComponent(0, lb - 2, 0, kb);
                comp_lk = comp_lk - lb2_fact * FactorialNN(0, 0, lb - 2, kb1) *
                                        t_tensor.at(t_cpn);
              }

              t_tensor.at(cpn) = comp_lk / FactorialNN(la, ka1, lb, kb1);

            }  //! recursive relation for lb >0 .and. la>0
            else {
              // initializing component
              double comp_lk{0};

              if (ka <= 2 * (la - 2)) {
                const int t_cpn = GetComponent(la - 2, lb, ka, kb);
                comp_lk = comp_lk + FactorialNN(la - 2, ka1, lb, kb1) *
                                        t_tensor.at(t_cpn);
              }

              if (kb <= 2 * (lb - 2)) {
                const auto l2_fact =
                    static_cast<double>((2.0 * la + lb - 1.0) / lb);
                const int t_cpn = GetComponent(la, lb - 2, ka, kb);
                comp_lk =
                    comp_lk - (l2_fact * FactorialNN(la, ka1, lb - 2, kb1)) *
                                  t_tensor.at(t_cpn);
              }

              for (int i = 1; i <= 3; i++) {
                int rk1;
                std::string rk2;
                tie(rk1, rk2) = NEta(coord.at(i - 1), kb1, kb2);
                const int rk_i = GetTensorComponent(lb, rk1, rk2);
                const double m = CoeffM(coord.at(i - 1), kb1, kb2);

                const auto l3_fact =
                    static_cast<double>((2.0 * (la + lb) - 1.0) / lb);
                const double const_fact =
                    l3_fact * m * FactorialNN(la, ka1, lb - 1, rk1);

                if (std::abs(const_fact) > DBL_EPSILON &&
                    rk_i <= 2 * (lb - 1)) {
                  const int t_cpn = GetComponent(la, lb - 1, ka, rk_i);
                  comp_lk =
                      comp_lk + const_fact * b.at(i - 1) * t_tensor.at(t_cpn);
                }
              }

              for (int i = 1; i <= 3; i++) {
                for (int j = 1; j <= 3; j++) {
                  const int n = 3 * (i - 1) + j;
                  const auto l4_fact =
                      static_cast<double>((2.0 * la - 1.0) / lb);
                  int rka1;
                  std::string rka2;
                  tie(rka1, rka2) = NEta(coord.at(i - 1), ka1, ka2);
                  int rkb1;
                  std::string rkb2;
                  tie(rkb1, rkb2) = NEta(coord.at(j - 1), kb1, kb2);

                  const int rk_i = GetTensorComponent(la, rka1, rka2);
                  const int rk_j = GetTensorComponent(lb, rkb1, rkb2);
                  const double m1 = CoeffM(coord.at(i - 1), ka1, ka2);
                  const double m2 = CoeffM(coord.at(j - 1), kb1, kb2);

                  const double const_factor =
                      l4_fact * m1 * m2 *
                      FactorialNN(la - 1, rka1, lb - 1, rkb1);
                  if (std::abs(const_factor) > DBL_EPSILON &&
                      rk_i <= 2 * (la - 1) && rk_j <= 2 * (lb - 1)) {
                    const int t_cpn = GetComponent(la - 1, lb - 1, rk_i, rk_j);
                    comp_lk = comp_lk +
                              const_factor * cc.at(n - 1) * t_tensor.at(t_cpn);
                  }
                }
              }

              t_tensor.at(cpn) = comp_lk / FactorialNN(la, ka1, lb, kb1);
            }
          }
          //*********************************************************
          cpn++;
        }
      }
    }
  }
  return t_tensor;
}