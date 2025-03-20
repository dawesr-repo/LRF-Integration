//
// Created by albpl on 3/20/2025.
//

#include "Tensor.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <tuple>


/// @brief ar,br,cab are the projection of the local molecular frame with
/// respect to the laboratory frame after rotation.
/// @return
std::vector<double> ar(const std::vector<double> &general_coordinates) {
  //{Az,Ax,Ay}
  const std::vector<double> &vec_ar = {
      cos(general_coordinates.at(1)),
      sin(general_coordinates.at(1)) * sin(general_coordinates.at(4)),
      cos(general_coordinates.at(4)) * sin(general_coordinates.at(1))};
  return vec_ar;
}

std::vector<double> br(const std::vector<double> &general_coordinates) {
  //{Bz,Bx,By}
  const std::vector<double> &vec_ar = {
      -cos(general_coordinates.at(2)),
      -sin(general_coordinates.at(2)) * sin(general_coordinates.at(5)),
      -cos(general_coordinates.at(5)) * sin(general_coordinates.at(2))};
  return vec_ar;
}

std::vector<double> cab(const std::vector<double> &general_coordinates) {
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

/// @param la, lb, ka, kb components of the T-tensor
/// @return linear component of the T-tensor
int getComponent(const int &la, const int &lb, const int &ka, const int &kb) {
  int cpn{0};
  for (int order = 1; order <= 15; ++order) {
    for (int lap = 0; lap <= order - 1; ++lap) {
      int lbp = order - lap - 1;
      for (int kap = 0; kap <= 2 * lap; ++kap) {
        for (int kbp = 0; kbp <= 2 * lbp; ++kbp) {
          if (la == lap && lb == lbp && ka == kap && kb == kbp) {
            return cpn;
          } else {
            cpn++;
          }
        }
      }
    }
  }
  return 0;
}



/// @brief Function that pass from linear index to the real-spherical index
/// @param i :linear index
/// @return component of the T-tensor in real-spherical notation
std::string getSplittingComponet(const int &i) {
  return i < 0 ? "-1" : i == 0 ? "0" : i % 2 == 1 ? "c" : "s";
}

/// @brief Function that pass from real-spherical index to the linear index
int getTensorComponent(const int &mult_ord, const int &k1,
                             const std::string &k2) {
  return (k1 < 0 || mult_ord < 0 || k1 > mult_ord) ? -1
         : k1 == 0                                 ? (k2 == "0" ? 0 : -1)
                   : (k2 == "s" ? 2 * k1 : (k2 == "c" ? 2 * k1 - 1 : 0));
}
/// @brief Auxiliar function for recursive relation function
std::tuple<int, std::string> nEta(const std::string &mu, const int &k1, const std::string &k2) {
  return mu == "x"   ? (k1 <= 1 ? std::make_tuple(0, "0") : make_tuple(k1 - 1, k2))
         : mu == "y" ? (k1 <= 1 ? std::make_tuple(0, "0")
                                : (k2 == "c" ? std::make_tuple(k1 - 1, "s")
                                             : std::make_tuple(k1 - 1, "c")))
                     : make_tuple(k1, k2);
}
/// @brief Auxiliar function for recursive relation function
double factorial(int n) {
  double res = 1.0;
  for (int i = 2; i <= n; i++) res *= static_cast<double>(i);
  return res;
}
/// @brief Auxiliar function for recursive relation function
double factorial_nn(const int &la, const int &ka1, const int &lb,
                          const int &kb1) {
  if (la < 0 || lb < 0 || ka1 < 0 || kb1 < 0 || ka1 > la || kb1 > lb) {
    return 0.0;
  } else {
    return sqrt((factorial(la + ka1) / factorial(la - ka1)) *
                (factorial(lb + kb1) / factorial(lb - kb1)));
  }
}
/// @brief Auxiliar function for recursive relation function
double coeffM(const std::string &mu, const int &k1, const std::string &k2) {
  double coeff_m = 0.0;

  if (mu == "x") {
    if (k1 == 1) {
      if (k2 == "c") {
        coeff_m = sqrt(2.0);
      }
    } else {
      coeff_m = static_cast<double>(k1);
    }
  } else if (mu == "y") {
    if (k1 == 1) {
      if (k2 == "s") {
        coeff_m = sqrt(2.0);
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



//******************* Tensor Methods *******************************************

/// @param system_dimension integer with the number of degrees of freedom of the system
/// @param coordinate_format format of the coordinates(supported formats:
/// "Euler_ZXZ", "Euler_ZYZ" and "Spherical")
/// @param coordinates R, beta1, beta2, alpha, gamma1, gamma2. R must be in
/// Angstroms and the angles in degrees.
Tensor:: Tensor(const int system_dimension,const std::vector <double> &coordinates, const std::string &coordinate_format) {

  system_dimension_ = system_dimension;
  coordinates_ = coordinates;
  coordinate_format_ = coordinate_format;
}

/// @brief Pass from user coordinates to general coordinates
/// @return the standard 6D vector in Euler_ZXZ vector.
///
std::vector<double> Tensor:: UserCoordinatesToGeneralCoordinates() const {
  // general_coodinates
  std::vector<double> vec(6, 0);
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
    std::cout << "Wrong dimimension: " << system_dimension_ << std::endl;
    //throw 1;
  }

  if (coordinate_format_ == "Euler_ZYZ") {
    vec.at(4) = vec.at(4) - 90.0;
    vec.at(5) = vec.at(5) - 90.0;
  } else if (coordinate_format_ == "Spherical") {
    vec.at(4) = 90.0 - vec.at(4);
    vec.at(5) = 90.0 - vec.at(5);
  }

  for (int i = 1; i < 6; ++i) {
    vec.at(i) = vec.at(i) * (M_PI / 180.0);
  }

  return vec;
}



/// @brief Calculate the T-tensor components for a given set of coordinates
/// @param max_t_tensor_order: Maximum order to calculate the T-Tensors. It must
///                            be a positive integer
///
/// @return 1D-vector with the T-tensor components (4-D tensor reshaped)

std::vector<double> Tensor::CalculateTensor(const int max_t_tensor_order) {

  std::vector<std::string> coord{"z", "x", "y"};  //! Cartesian Axis Labels

  if (coordinates_.size() != system_dimension_) {
    std::cout << "coordinates size must be equal to dimension" << std::endl;
    //throw 0;
  }
  if (coordinates_.at(0) < 1) {
    std::cout << "coordinates size must be equal to dimension" << std::endl;
    //throw 0;
  }

  const auto general_coordinates_zxz =  UserCoordinatesToGeneralCoordinates();

  auto a = br(general_coordinates_zxz);
  auto b = br(general_coordinates_zxz);
  auto cc = cab(general_coordinates_zxz);

  const int cpns_per_order[15]{1,    7,    26,   70,   155,  301,  532, 876,
                               1365, 2035, 2926, 4082, 5551, 7385, 9640};

  std::vector<double> t_tensor(cpns_per_order[max_t_tensor_order - 1], 0);

  int cpn{0};

  for (int order = 1; order <= max_t_tensor_order; ++order) {
    for (int la = 0; la <= order - 1; ++la) {
      int lb = order - la - 1;
      for (int ka = 0; ka <= 2 * la; ++ka) {
        for (int kb = 0; kb <= 2 * lb; ++kb) {
          // Calculating T-Tensor Component
          const std::string ka2 = getSplittingComponet(ka);
          const std::string kb2 = getSplittingComponet(kb);

          const int ka1 = floor((ka + 1.0) / 2.0);
          const int kb1 = floor((kb + 1.0) / 2.0);

          if (la == 0 && lb == 0) {  // 1st order tensor component
            t_tensor.at(0) = 1.0;
          } else {
            // recursive relation for lb = 0
            if (lb == 0) {
              // initializating component
              double comp_lk{0};
              const auto la_fact = static_cast<double>((2.0 * la - 1.0) / la);

              // loop though every coodinate axis
              for (int i = 1; i <= 3; i++) {
                // new multipole components
                int rk1;
                std::string rk2;
                tie(rk1, rk2) = nEta(coord.at(i - 1), ka1, ka2);
                const int rk_ = getTensorComponent(la, rk1, rk2);
                const double m = coeffM(coord.at(i - 1), ka1, ka2);
                // coeffcient NN of the recurence
                const double fact_nn = factorial_nn(la - 1, rk1, 0, 0);

                if (std::abs(m) > DBL_EPSILON && la >= 1 && fact_nn > DBL_EPSILON &&
                    rk_ <= 2 * (la - 1)) {
                  const int t_cpn = getComponent(la - 1, 0, rk_, 0);
                  double prod_comp = a.at(i - 1) * t_tensor.at(t_cpn);
                  double fact_prod = la_fact * m * fact_nn;
                  comp_lk = comp_lk + fact_prod * prod_comp;
                }
              }
              // Second term of the recurence
              if (la >= 2 && ka <= 2 * (la - 2) && ka >= 0) {
                const auto la2_fact = static_cast<double>((la - 1.0) / la);

                const int t_cpn = getComponent(la - 2, 0, ka, 0);
                comp_lk = comp_lk - la2_fact * factorial_nn(la - 2, ka1, 0, 0) *
                                        t_tensor.at(t_cpn);
              }

              t_tensor.at(cpn) = comp_lk / factorial_nn(la, ka1, lb, kb1);
            }
            // recursive relation for la = 0
            else if (la == 0) {
              // initializating component
              double comp_lk{0};
              const auto lb_fact = static_cast<double>((2.0 * lb - 1.0) / lb);

              // loop though every coodinate axis
              for (int i = 1; i <= 3; i++) {
                // new multipole components
                int rk1;
                std::string rk2;
                tie(rk1, rk2) = nEta(coord.at(i - 1), kb1, kb2);
                const int rk_ = getTensorComponent(lb, rk1, rk2);
                const double m = coeffM(coord.at(i - 1), kb1, kb2);
                // coeffcient NN of the recurence
                const double fact_nn = factorial_nn(0, 0, lb - 1, rk1);

                if (std::abs(m) > DBL_EPSILON && lb >= 1 && fact_nn > DBL_EPSILON &&
                    rk_ <= 2 * (lb - 1) && rk_ >= 0) {
                  const int t_cpn = getComponent(0, lb - 1, 0, rk_);
                  comp_lk = comp_lk + lb_fact * m * fact_nn * b.at(i - 1) *
                                          t_tensor.at(t_cpn);
                }
              }
              // Second term of the recurence
              if (lb >= 2 && kb <= 2 * (lb - 2) && kb >= 0) {
                const auto lb2_fact = static_cast<double>((lb - 1.0) / lb);
                const int t_cpn = getComponent(0, lb - 2, 0, kb);
                comp_lk = comp_lk - lb2_fact * factorial_nn(0, 0, lb - 2, kb1) *
                                        t_tensor.at(t_cpn);
              }

              t_tensor.at(cpn) = comp_lk / factorial_nn(la, ka1, lb, kb1);

            }  //! recursive relation for lb >0 .and. la>0
            else {
              // initializating component
              double comp_lk{0};

              if (ka <= 2 * (la - 2)) {
                const int t_cpn = getComponent(la - 2, lb, ka, kb);
                comp_lk = comp_lk + factorial_nn(la - 2, ka1, lb, kb1) *
                                        t_tensor.at(t_cpn);
              }

              if (kb <= 2 * (lb - 2)) {
                const auto l2_fact =
                    static_cast<double>((2.0 * la + lb - 1.0) / lb);
                const int t_cpn = getComponent(la, lb - 2, ka, kb);
                comp_lk =
                    comp_lk - (l2_fact * factorial_nn(la, ka1, lb - 2, kb1)) *
                                  t_tensor.at(t_cpn);
              }

              for (int i = 1; i <= 3; i++) {
                int rk1;
                std::string rk2;
                tie(rk1, rk2) = nEta(coord.at(i - 1), kb1, kb2);
                const int rk_i = getTensorComponent(lb, rk1, rk2);
                const double m = coeffM(coord.at(i - 1), kb1, kb2);

                const auto l3_fact =
                    static_cast<double>((2.0 * (la + lb) - 1.0) / lb);
                const double const_fact =
                    l3_fact * m * factorial_nn(la, ka1, lb - 1, rk1);

                if (std::abs(const_fact) > DBL_EPSILON && rk_i <= 2 * (lb - 1)) {
                  const int t_cpn = getComponent(la, lb - 1, ka, rk_i);
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
                  tie(rka1, rka2) = nEta(coord.at(i - 1), ka1, ka2);
                  int rkb1;
                  std::string rkb2;
                  tie(rkb1, rkb2) = nEta(coord.at(j - 1), kb1, kb2);

                  const int rk_i = getTensorComponent(la, rka1, rka2);
                  const int rk_j = getTensorComponent(lb, rkb1, rkb2);
                  const double m1 = coeffM(coord.at(i - 1), ka1, ka2);
                  const double m2 = coeffM(coord.at(j - 1), kb1, kb2);

                  const double const_factor =
                      l4_fact * m1 * m2 *
                      factorial_nn(la - 1, rka1, lb - 1, rkb1);
                  if (std::abs(const_factor) > DBL_EPSILON && rk_i <= 2 * (la - 1) &&
                      rk_j <= 2 * (lb - 1)) {
                    const int t_cpn = getComponent(la - 1, lb - 1, rk_i, rk_j);
                    comp_lk = comp_lk +
                              const_factor * cc.at(n - 1) * t_tensor.at(t_cpn);
                  }
                }
              }

              t_tensor.at(cpn) = comp_lk / factorial_nn(la, ka1, lb, kb1);
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