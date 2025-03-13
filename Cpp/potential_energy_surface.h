//
// Created by albpl on 3/12/2025.
//

#ifndef PES_H
#define PES_H

#include <iostream>
#include <utility>
#include <vector>
using namespace std;


class PotentialEnergySurface {
  string file_name_;

  const double C1 = 627.5095;
  const double C2 = 0.529177249;
  const double C3 = 349.757;

  int max_t_tensor_order_;
  int m_fit_[15];
  int i_fit_[15];
  int d_fit_[15];
  // Zero!
  double zero_;

  // Multipoles!
  double a_mult_[225] ;
  double b_mult_[225] ;

  // Polarizability!
  vector<double> a_pol_[6][12];
  vector<double> b_pol_[6][12];

  //
  // Dispersion
  vector<double> disp_[5][10][5][10];


public:
  explicit PotentialEnergySurface(const string &filename);
  const double GetAsymptote()const {return zero_;};
  const double EvaluateLRF( const int &dim,
                      const vector <double> &coordinates,
                      const string &coordinate_format);
private:
  const double MultipoleInteraction(const double &r,const vector <double> &t_tensors);
  const double MultipoleOrder(const int &order,const vector <double> &t_tensors);

  const double InductionInteraction(const double &r,const vector <double> &t_tensors);
  const double InductionOrder(const int &order,const double &r,
                              const vector <double> &t_tensors,const int &index);
  const double InductionComponent( const int &i,const int &j,
                                    const int &l1,const int &l2,
                                    const vector <double> &t_tensors,
                                    const int &index);

  const double DispersionInteraction(const double &r,const vector <double> &t_tensors);
  const double DispersionOrder(const double &r,
                               const vector <double> &t_tensors,
                               const int &order);
  const double DispersionComponent(const int &l1,const int &l2,
                                   const int &t1,const int &t2,
                                   const vector <double> &t_tensors);
  const vector <double> CalculateTensor(const vector <double>  &general_coordinates);
  const double GetTotalInteractionEnergy(const vector <double>  &general_coordinates);
  bool ReadParameters();
};

#endif //PES_H
