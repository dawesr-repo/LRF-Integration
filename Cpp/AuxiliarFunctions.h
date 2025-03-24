//
// Created by albpl on 3/24/2025.
//

#ifndef AUXILIAR_FUNCTIONS_H
#define AUXILIAR_FUNCTIONS_H
#include <fstream>
#include <limits>

template<typename T>
class AuxiliarFunctions {
  /// @brief function that skips n lines in the file fp
  /// @param fp   ifstream object
  /// @param n_lines integer with the number of lines to skip
public:
  static void skipLines(std::ifstream &fp,int n_lines);
  static void readArray(std::ifstream &fp, T *arr,int n);
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

#endif //AUXILIAR_FUNCTIONS_H
