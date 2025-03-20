//
// Created by albpl on 3/19/2025.
//

#include "PotentialEnergySurface.h"
#include "Tensor.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <limits>

//************ Auxiliar Functions Definition **********************

/// @brief function that skips nlines in the file fp
/// @param fp   ifstream object
/// @param nlines integer with the number of lines to skip
void skipLines(std::ifstream &fp, int nlines) {
  for (int i = 0; i < nlines; ++i) {
    fp.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
}

/// @brief function that reads an array of n elements from the file fp
/// @param fp   ifstream object
/// @param arr  pointer to the array where the elements will be stored
/// @param n    integer with the number of elements to read
template <class T>
void readArray(std::ifstream &fp, T *arr, int n) {
  for (int i = 0; i < n; ++i) {
    fp >> arr[i];
  }
}


//************ PotentialEnergySurface Methods Definition **********************

/// @brief Constructor of the class PotentialEnergySurface
/// @param filename: string with the name of the file that contains the
/// parameters of the potential energy surface.
PotentialEnergySurface::PotentialEnergySurface(const std::string &filename)
     : file_name_(filename) {
  bool err = ReadParameters();

  if (!err) {
    std::cerr << "Failed to ReadParameters from " << file_name_ << std::endl;
  }
}
/// @brief ReadParameters reads the multipole, polarizability and dispersion
/// coefficients from the file exported
///        from MATLAB and stores them in the class PotentialEnergySurface. Also
///        it calculates the maximum order of the t-tensors that will be used in
///        the calculations.
/// @return bool: true if the parameters were read successfully, false
/// otherwise.
bool PotentialEnergySurface::ReadParameters() {
  std::ifstream infile(file_name_);
  std::string singleline, tmp;
  if (!infile) {
    std::cerr << "Failed to open file: " << file_name_ << std::endl;
    return false;
  }

  skipLines(infile, 8);

  infile >> tmp >> singleline;
  char* end = nullptr;
  zero_ = strtod(singleline.c_str(),&end);

  skipLines(infile, 3);
  // Multipoles
  infile >> tmp;
  readArray(infile, m_fit_, 15);
  infile >> tmp;
  readArray(infile, i_fit_, 15);
  infile >> tmp;
  readArray(infile, d_fit_, 15);
  infile >> tmp;
  readArray(infile, a_mult_, 225);
  infile >> tmp;
  readArray(infile, b_mult_, 225);

  const int iord = *std::max_element(i_fit_, i_fit_ + 15);
  const int temp = *std::max_element(m_fit_, m_fit_ + 15);
  const int mord = std::max(temp, iord - 3);
  const int dord = *std::max_element(d_fit_, d_fit_ + 15);

  max_t_tensor_order_ = std::max(std::max(iord - 2, mord), dord - 3);

  if (iord >= 4) {
    for (int i = 1; i <= iord - 3; ++i) {
      for (int j = i; j <= iord - 3; ++j) {
        if (i + j <= iord - 2) {
          const int ln = (2 * i + 1) * (2 * j + 1);
          std::vector<double> &arr1 = a_pol_[i - 1][j - 1];
          arr1.resize(ln, 0);
          infile >> tmp;
          readArray(infile, arr1.data(), ln);

          std::vector<double> &arr2 = b_pol_[i - 1][j - 1];
          arr2.resize(ln, 0);
          infile >> tmp;
          readArray(infile, arr2.data(), ln);
        }
      }
    }
  }

  if (dord >= 6) {
    for (int l1 = 1; l1 <= dord - 5; ++l1) {
      for (int l2 = l1; l2 <= dord - 5; ++l2) {
        for (int t1 = 1; t1 <= dord - 5; ++t1) {
          for (int t2 = t1; t2 <= dord - 5; ++t2) {
            if (l1 + l2 + t1 + t2 <= dord - 2) {
              const int ln =
                  (2 * l1 + 1) * (2 * l2 + 1) * (2 * t1 + 1) * (2 * t2 + 1);
              std::vector<double> &arr = disp_[l1 - 1][l2 - 1][t1 - 1][t2 - 1];
              arr.resize(ln, 0);
              infile >> tmp;
              readArray(infile, arr.data(), ln);
            }
          }
        }
      }
    }
  }

  infile.close();
  return true;
}


//************ PotentialEnergySurface Methods Definition **********************

/// @brief EvaluateLRF calculates the energy of a system with a given the given
/// set of coordinates by the user
/// @param dim integer with the number of degrees of freedom of the system
/// @param coordinates vector with the coordinates of the system
/// @param coordinate_format format of the coordinates(supported formats:
/// "Euler_ZXZ", "Euler_ZYZ" and "Spherical")
/// @return double with the total interaction energy
double PotentialEnergySurface::EvaluateLRF(
    const int &dim, const std::vector<double> &coordinates,
    const std::string &coordinate_format) {


  auto* t_tensor = new Tensor(dim,coordinates,coordinate_format);
  const std::vector<double> t = t_tensor->CalculateTensor(max_t_tensor_order_);
  const double r = t_tensor->GetIntermolecularDistance();

  return MultipoleInteraction(r, t) +
         InductionInteraction(r, t) +
         DispersionInteraction(r, t);


}