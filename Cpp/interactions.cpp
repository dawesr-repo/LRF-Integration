//
// Created by albpl on 3/13/2025.
//
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <tuple>

#include "potential_energy_surface.h"
/// @brief component of the T-tensor
/// @param la,lb,ka,kb components of the T-tensor
/// @return Linear index component of the T-tensor
int getComponents(const int &la, const int &lb, const int &ka, const int &kb) {
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

//************ PotentialEnergySurface Methods Definition **********************

//********************  MULTIPOLE INTERACTION  ********************************
/// @brief Calculate the multipole interaction between two molecules by order
/// @param r distance between the two molecules
/// @param t_tensors vector with the T-tensors components
/// @return the multipole interaction between the two molecules
const double PotentialEnergySurface::MultipoleOrder(
    const int &order, const vector<double> &t_tensors) {
  double multipole_order{0.0};
  const double EPS = DBL_EPSILON;

  for (int i = 0; i <= order - 1; ++i) {
    const int j = order - 1 - i;

    for (int ci = 0; ci <= 2 * i; ++ci) {
      const double Qai = a_mult_[i * i + ci];

      if (abs(Qai) > EPS) {
        for (int cj = 0; cj <= 2 * j; ++cj) {
          const double Qbj = a_mult_[j * j + cj];

          if (abs(Qbj) > EPS) {
            const int t_cpn = getComponents(i, j, ci, cj);
            multipole_order = multipole_order + Qai * Qbj * t_tensors.at(t_cpn);
          }
        }
      }
    }
  }
  return multipole_order;
}
/// @brief calculate the multipole interaction between two molecules
/// @param r distance between the two molecules
/// @param t_tensors vector with the T-tensors components
/// @return the multipole interaction between the two molecules
const double PotentialEnergySurface::MultipoleInteraction(
    const double &r, const vector<double> &t_tensors) {
  double multipole_sph = 0.0;
  for (int order = 1; order <= 15; ++order) {
    if (m_fit_[order - 1] > 0) {
      multipole_sph = multipole_sph + (C3 * C1 * pow(C2, order)) *
                                          MultipoleOrder(order, t_tensors) /
                                          pow(r, order);
    }
  }

  return multipole_sph;
}

//********************  INDUCTION INTERACTION  ********************************

/// @brief Calculate the induction interaction between two molecules by order
/// @param r distance between the two molecules
/// @param t_tensors vector with the T-tensors components
/// @return the induction interaction between the two molecules

const double PotentialEnergySurface::InductionInteraction(
    const double &r, const vector<double> &t_tensors) {
  double induction_sph = 0.0;
  for (int order = 1; order <= 15; ++order) {
    if (i_fit_[order - 1] > 0) {
      induction_sph =
          induction_sph +
          InductionOrder(order, r, t_tensors, 1);  // induction of B over A
      +InductionOrder(order, r, t_tensors, 0);     // induction of A over B
    }
  }

  return induction_sph;
}

/// @brief Calculate the induction interaction between two molecules by order
/// @param order order of the induction interaction
/// @param r distance between the two molecules
/// @param t_tensors vector with the T-tensors components
/// @param index indicate if Im calculating pol over A or pol over B
/// @return the induction interaction between the two molecules

const double PotentialEnergySurface ::InductionOrder(
    const int &order, const double &r, const vector<double> &t_tensors,
    const int &index) {
  double res{0.0};

  for (int l1 = 1; l1 <= order - 3; ++l1) {
    for (int l2 = 1; l2 <= order - 3; ++l2) {
      if (l1 + l2 + 2 <= order) {
        for (int i = 0; i <= order - 2 - l1 - l2; ++i) {
          for (int j = 0; j <= order - 2 - l1 - l2; ++j) {
            if (i + j + l1 + l2 + 2 == order) {
              res = res + InductionComponent(i, j, l1, l2, t_tensors, index);
            }
          }
        }
      }
    }
  }

  return (-0.5 * (C3 * C1 * pow(C2, order)) * res) / (pow(r, order));
}

/// @brief Calculate the induction interaction between two molecules by
/// components
/// @param i,j,l1,l2 components of the induction interaction
/// @param t_tensors vector with the T-tensors components
/// @param index indicate if Im calculating pol over A or pol over B
/// @return the induction interaction between the two molecules

const double PotentialEnergySurface ::InductionComponent(
    const int &i, const int &j, const int &l1, const int &l2,
    const vector<double> &t_tensors, const int &index) {
  double res{0.0};
  const int ni = 2 * i + 1;
  const int nj = 2 * j + 1;
  const int nl1 = 2 * l1 + 1;
  const int nl2 = 2 * l2 + 1;
  const int lmin = min(l1, l2);
  const int lmax = max(l1, l2);
  const double EPS = DBL_EPSILON;

  const auto mult_cpn = index == 1 ? a_mult_ : b_mult_;
  const auto pol_arr =
      index == 1 ? b_pol_[lmin - 1][lmax - 1] : a_pol_[lmin - 1][lmax - 1];

  for (int ci = 1; ci <= ni; ++ci) {
    const double qai = mult_cpn[i * i + ci - 1];
    if (abs(qai) > EPS) {
      for (int cj = 1; cj <= nj; ++cj) {
        const double qbj = mult_cpn[j * j + cj - 1];
        if (abs(qbj) > EPS) {
          for (int k1 = 1; k1 <= nl1; ++k1) {
            for (int k2 = 1; k2 <= nl2; ++k2) {
              const int cpn = l1 > l2 ? (k2 - 1) * (2 * l1 + 1) + k1
                                      : (k1 - 1) * (2 * l2 + 1) + k2;
              const double comp_a_k1_k2 = pol_arr.at(cpn - 1);

              if (abs(comp_a_k1_k2) > EPS) {
                // index indicate if Im calculating pol over A
                //  or pol over B
                if (index == 0) {
                  const int t_cpn_1 = getComponents(l1, i, k1 - 1, ci - 1);
                  const int t_cpn_2 = getComponents(l2, j, k2 - 1, cj - 1);
                  res =
                      res + qai * qbj * comp_a_k1_k2 *
                                (t_tensors.at(t_cpn_1) * t_tensors.at(t_cpn_2));

                } else {
                  const int t_cpn_1 = getComponents(i, l1, ci - 1, k1 - 1);
                  const int t_cpn_2 = getComponents(j, l2, cj - 1, k2 - 1);
                  res =
                      res + qai * qbj * comp_a_k1_k2 *
                                (t_tensors.at(t_cpn_1) * t_tensors.at(t_cpn_2));
                }
              }
            }
          }
        }
      }
    }
  }

  return res;
}

// *********************  DISPERSION INTERACTION  ****************************

/// @brief Calculate the dispersion interaction between two molecules by order
/// @param r distance between the two molecules
/// @param t_tensors vector with the T-tensors components
/// @return the dispersion interaction between the two molecules

const double PotentialEnergySurface::DispersionInteraction(
    const double &r, const vector<double> &t_tensors) {
  double dispersion_sph{0.0};

  for (int order = 1; order <= 15; ++order) {
    if (d_fit_[order - 1] > 0) {
      dispersion_sph = dispersion_sph + DispersionOrder(r, t_tensors, order);
    }
  }
  return dispersion_sph;
}

const double PotentialEnergySurface::DispersionOrder(
    const double &r, const vector<double> &t_tensors, const int &order) {
  double res{0.0};

  for (int l1 = 1; l1 <= order - 2; ++l1) {
    for (int l2 = 1; l2 <= order - 2 - l1; ++l2) {
      for (int t1 = 1; t1 <= order - 2 - l1 - l2; ++t1) {
        for (int t2 = 1; t2 <= order - 2 - l1 - l2 - t1; ++t2) {
          if (l1 + l2 + t1 + t2 + 2 == order) {
            res = res + DispersionComponent(l1, l2, t1, t2, t_tensors);
          }
        }
      }
    }
  }

  return -((C3 * C1 * pow(C2, order)) * res) / pow(r, order);
}

const double PotentialEnergySurface::DispersionComponent(
    const int &l1, const int &l2, const int &t1, const int &t2,
    const vector<double> &t_tensors) {
  double res{0.0};
  const double EPS = DBL_EPSILON;
  const int lmin = min(l1, l2);
  const int lmax = max(l1, l2);
  const int tmin = min(t1, t2);
  const int tmax = max(t1, t2);
  const auto disp_arr = disp_[lmin - 1][lmax - 1][tmin - 1][tmax - 1];

  for (int li = 0; li <= 2 * l1; ++li) {
    for (int lj = 0; lj <= 2 * l2; ++lj) {
      for (int ti = 0; ti <= 2 * t1; ++ti) {
        for (int tj = 0; tj <= 2 * t2; ++tj) {
          const int cpn =
              (l1 > l2 && t1 > t2)
                  ? lj * (2 * l1 + 1) * (2 * t2 + 1) * (2 * t1 + 1) +
                        li * (2 * t2 + 1) * (2 * t1 + 1) + tj * (2 * t1 + 1) +
                        ti + 1
              : (l1 > l2 && t1 <= t2)
                  ? lj * (2 * l1 + 1) * (2 * t1 + 1) * (2 * t2 + 1) +
                        li * (2 * t1 + 1) * (2 * t2 + 1) + ti * (2 * t2 + 1) +
                        tj + 1
              : (l1 <= l2 && t1 > t2)
                  ? li * (2 * l2 + 1) * (2 * t2 + 1) * (2 * t1 + 1) +
                        lj * (2 * t2 + 1) * (2 * t1 + 1) + tj * (2 * t1 + 1) +
                        ti + 1
                  : li * (2 * l2 + 1) * (2 * t1 + 1) * (2 * t2 + 1) +
                        lj * (2 * t1 + 1) * (2 * t2 + 1) + ti * (2 * t2 + 1) +
                        tj + 1;

          const double disp_coeff = disp_arr.at(cpn - 1);

          if (abs(disp_coeff) > EPS) {
            const int t_cpn_1 = getComponents(l1, t1, li, ti);
            const int t_cpn_2 = getComponents(l2, t2, lj, tj);

            res = res +
                  disp_coeff * t_tensors.at(t_cpn_1) * t_tensors.at(t_cpn_2);
          }
        }
      }
    }
  }

  return res;
}