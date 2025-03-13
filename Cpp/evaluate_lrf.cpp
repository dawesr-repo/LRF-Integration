//
// Created by albpl on 3/12/2025.
//


#include <fstream>
#include "potential_energy_surface.h"
#include <limits>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cfloat>



const vector <double> userCoordinatesToGeneralCoordinates(const int &dim,
                                                          const string &coordinate_format,
                                                          const vector <double>  &user_coordinates){
  //general_coodinates
  vector<double> vec(6,0) ;
  if (dim>=2 && dim<=6){

    vec.at(0) = user_coordinates.at(0);  //R
    vec.at(1) = user_coordinates.at(1);  //beta1
    if (dim==3){
      vec.at(4) = user_coordinates.at(2);//gamma1
    }else if (dim>=4){
      vec.at(2) = user_coordinates.at(2);//beta2
      vec.at(3) = user_coordinates.at(3);//alpha
    }
    if (dim==5){
      vec.at(4) = user_coordinates.at(4);//gamma1
    }
    if (dim==6){
      vec.at(5) = user_coordinates.at(5);//gamma2
    }
  }else{
    cout<<"Wrong dimimension: "<<dim<<endl;
    throw 0;
  }

  if (coordinate_format == "Euler_ZYZ"){
    vec.at(4) = vec.at(4) - 90.0;
    vec.at(5) = vec.at(5) - 90.0;
  }else if (coordinate_format == "Spherical"){
    vec.at(4) = 90.0 - vec.at(4);
    vec.at(5) = 90.0 - vec.at(5);
  }
  double pi = M_PI;
  for (int i = 1; i < 6; ++i) {
    vec.at(i) = vec.at(i)*(M_PI / 180.0);
  }

  const vector <double> &general_coordinates_zxz = {vec};
  return general_coordinates_zxz;
}


const double PotentialEnergySurface::GetTotalInteractionEnergy(const vector <double>  &general_coordinates){

  const vector <double>  t_tensors = CalculateTensor(general_coordinates);
  const double r = general_coordinates.at(0);
  return  MultipoleInteraction(r,t_tensors) +
          InductionInteraction(r,t_tensors) +
          DispersionInteraction(r,t_tensors);
}





//************ PotentialEnergySurface Methods Definition **********************

const double PotentialEnergySurface::EvaluateLRF( const int &dim,
                      const vector <double> &coordinates,
                      const string &coordinate_format){

  if (coordinates.size() != dim){
    cout<<"coordinates size must be equal to dimimension"<<endl;
    throw 0;
  }
  if (coordinates.at(0)<1) {
    cout<<"coordinates size must be equal to dimimension"<<endl;
    throw 0;
  }

  auto general_coordinates_zxz = userCoordinatesToGeneralCoordinates(dim,
                                                                     coordinate_format,
                                                                     coordinates);

  return GetTotalInteractionEnergy(general_coordinates_zxz);
}
