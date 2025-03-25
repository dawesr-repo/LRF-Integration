//
// Created by albpl on 3/24/2025.
//

#ifndef AUXILIAR_FUNCTIONS_H
#define AUXILIAR_FUNCTIONS_H
#include <fstream>
#include <limits>
#include <cmath>
#include <numeric>


template<typename T>
class AuxiliarFunctions {
  /// @brief function that skips n lines in the file fp
  /// @param fp   ifstream object
  /// @param n_lines integer with the number of lines to skip
public:
  static void skipLines(std::ifstream &fp,int n_lines);
  static void readArray(std::ifstream &fp, T *arr,int n);
  static bool areEqual(const T *arr1,const T *arr2, int n);
  static int getComponents(const int &la, const int &lb, const int &ka, const int &kb);
  static int getComponents_v2(const int &la, const int &lb, const int &ka, const int &kb);
};

template<typename T>
/// @brief function that skips n lines in the file fp
/// @param fp   ifstream object
/// @param n_lines integer with the number of lines to skip
void AuxiliarFunctions<T>::skipLines(std::ifstream &fp, const int n_lines) {
  for (int i = 0; i < n_lines; ++i) {
    fp.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
}

/// @brief function that reads an array of n elements from the file fp
/// @param fp   ifstream object
/// @param arr  pointer to the array where the elements will be stored
/// @param n    integer with the number of elements to read
template <class T>
void AuxiliarFunctions<T>::readArray(std::ifstream &fp, T *arr, const int n) {

  for (int i = 0; i < n; ++i) {
    fp >> arr[i];
  }
}


/// @brief function that reads an array of n elements from the file fp
/// @param arr1   array of datatype T
/// @param arr2   array of datatype T
/// @param n    integer with the number of elements in the array.
/// @return true or false whether the values of arr1 and arr2 are equal.
template <class T>
bool AuxiliarFunctions<T>::areEqual(const T *arr1,const T *arr2, const int n) {

  const double eps = pow(10,-8);
  for (int i = 0; i < n; ++i) {
    if (std::abs(arr1[i]-arr2[i])>eps)
    return false;
  }
  return true;
}

/// @brief function that reads an array of n elements from the file fp
/// @param la   T_tensor component
/// @param ka   T_tensor component
/// @param lb    integer with the number of elements in the array.
/// @param kb   T_tensor component
/// @return the equivalent component in a linear array.
template <class T>
int AuxiliarFunctions<T>::getComponents(const int &la, const int &lb, const int &ka, const int &kb) {
  int cpn{0};
  for (int order = 1; order <= 15; ++order) {
    for (int lap = 0; lap <= order - 1; ++lap) {
      int lbp = order - lap - 1;
      for (int kap = 0; kap <= 2 * lap; ++kap) {
        for (int kbp = 0; kbp <= 2 * lbp; ++kbp) {
          if (la == lap && lb == lbp && ka == kap && kb == kbp) {
            return cpn;
          }
          cpn++;

        }
      }
    }
  }
  return 0;
}
/// @brief function that reads an array of n elements from the file fp
/// @param la   T_tensor component
/// @param ka   T_tensor component
/// @param lb    integer with the number of elements in the array.
/// @param kb   T_tensor component
/// @return the equivalent component in a linear array.
template <class T>
int AuxiliarFunctions<T>::getComponents_v2(const int &la, const int &lb, const int &ka, const int &kb) {
  const int cpns_per_order[16]{0, 1,    7,    26,   70,   155,  301,  532, 876,
                               1365, 2035, 2926, 4082, 5551, 7385, 9640};
  const int order = la+lb+1;
  return  cpns_per_order[order-1] +
            2*(order-1)*(la-1)*la +
            (2*order-1)*la -
            2*(la-1)*la*(2*la-1)/3 +
           ka*(2*lb+1) + kb  ;

}


#endif //AUXILIAR_FUNCTIONS_H
