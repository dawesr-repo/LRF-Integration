//
// Created by albpl on 3/19/2025.
//

#ifndef POTENTIAL_ENERGY_SURFACE_H
#define POTENTIAL_ENERGY_SURFACE_H

#include <iostream>
#include <vector>

class PotentialEnergySurface {
  std::string file_name_;

  const double C1_ = 627.5095;
  const double C2_ = 0.529177249;
  const double C3_ = 349.755088236337;

  int max_t_tensor_order_;
  int m_fit_[15];
  int i_fit_[15];
  int d_fit_[15];
  // Zero!
  double zero_;

  // Multipoles!
  double a_mult_[225];
  double b_mult_[225];

  // Polarizability!
  std::vector<double> a_pol_[6][12];
  std::vector<double> b_pol_[6][12];

  //
  // Dispersion
  std::vector<double> disp_[5][10][5][10];

 public:
  explicit PotentialEnergySurface(const std::string &filename);
  [[nodiscard]] double GetAsymptote() const { return zero_; };
  [[nodiscard]] double EvaluateLRF(const int &dim, const std::vector<double> &coordinates,
                     const std::string &coordinate_format) const;
  [[nodiscard]] double MultipoleInteraction(const double &r,
                            const std::vector<double> &t_tensors) const;
  [[nodiscard]] double InductionInteraction(const double &r,
                              const std::vector<double> &t_tensors) const;
  [[nodiscard]] double DispersionInteraction(const double &r,
                             const std::vector<double> &t_tensors) const;

  [[nodiscard]] const double *getMultipoleCoefficients(const std::string &label) const {return label=="A"? a_mult_:b_mult_;}
 private:

  [[nodiscard]] double MultipoleOrder(const int &order, const std::vector<double> &t_tensors) const;

  [[nodiscard]] double InductionOrder(const int &order, const double &r,
                        const std::vector<double> &t_tensors, const int &index) const;
  [[nodiscard]] double InductionComponent(const int &i, const int &j, const int &l1,
                            const int &l2, const std::vector<double> &t_tensors,
                            const int &index) const;
  [[nodiscard]] double DispersionOrder(const double &r, const std::vector<double> &t_tensors,
                         const int &order) const;
  [[nodiscard]] double DispersionComponent(const int &l1, const int &l2, const int &t1,
                             const int &t2,
                             const std::vector<double> &t_tensors) const;

  bool ReadParameters();

};

#endif  // POTENTIAL_ENERGY_SURFACE_H
