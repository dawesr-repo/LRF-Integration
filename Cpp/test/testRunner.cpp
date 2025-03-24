
//
// Created by albpl on 3/19/2025.
//

#include <gtest/gtest.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "../AuxiliarFunctions.h"
#include "../Tensor.h"

constexpr int system_dimension = 6;
const std::vector <double> coordinates {10.00,10.0,20.0,30.0,40.0,50.0};
const std::string coordinate_format = "Euler_ZYZ";

auto* t_tensor = new Tensor(system_dimension,coordinates,coordinate_format);

TEST(Tensor, Factorial) {

  EXPECT_EQ(0, t_tensor->Factorial(-1));
  EXPECT_EQ(1, t_tensor->Factorial(0));
  EXPECT_EQ(1, t_tensor->Factorial(1));
  EXPECT_EQ(24, t_tensor->Factorial(4));

}

TEST(Tensor, FactorialNN) {

  EXPECT_EQ(0, t_tensor->FactorialNN(-1,1,1,1));
  EXPECT_EQ(0, t_tensor->FactorialNN(1,-1,1,1));
  EXPECT_EQ(0, t_tensor->FactorialNN(1,1,-1,1));
  EXPECT_EQ(0, t_tensor->FactorialNN(1,1,1,-1));
  EXPECT_EQ(0, t_tensor->FactorialNN(1,2,1,1));
  EXPECT_EQ(0, t_tensor->FactorialNN(1,1,1,2));
  EXPECT_EQ(1, t_tensor->FactorialNN(0,0,0,0));
  EXPECT_EQ(2, t_tensor->FactorialNN(1,1,1,1));

}

TEST(Tensor, GetSplittingComponent) {

  EXPECT_EQ("-1", t_tensor->GetSplittingComponent(-1));
  EXPECT_EQ("0", t_tensor->GetSplittingComponent(0));
  EXPECT_EQ("c", t_tensor->GetSplittingComponent(1));
  EXPECT_EQ("s", t_tensor->GetSplittingComponent(2));

}

TEST(Tensor, GetTensorComponent) {

  EXPECT_EQ(-1, t_tensor->GetTensorComponent(-1,0,"0"));
  EXPECT_EQ(-1, t_tensor->GetTensorComponent(1,2,"0"));
  EXPECT_EQ(-1, t_tensor->GetTensorComponent(0,0,"c"));
  EXPECT_EQ(0, t_tensor->GetTensorComponent(4,0,"0"));
  EXPECT_EQ(1, t_tensor->GetTensorComponent(4,1,"c"));
  EXPECT_EQ(4, t_tensor->GetTensorComponent(4,2,"s"));

}

TEST(Tensor, Ar_Components) {

  std::vector<double> from_text(9,0);
  std::ifstream infile("../../../testing_datafiles/t_tensors/Ar.txt",std::ios_base::in);


  if (!infile) {
    std::cerr << " * Failed to open file "  << std::endl;
    if (infile.fail()) {
      // Print a more detailed error message using
      // strerror
      std::cerr << "Error details: " << strerror(errno)
           << std::endl;
    }
  }

  AuxiliarFunctions<double>::readArray(infile, from_text.data(), 9);

  const auto user_coordinates =  std::vector<double>(from_text.begin(),from_text.begin()+6);
  const auto ar_test =
      std::vector<double>(from_text.begin() + 6, from_text.end());

  const auto *t_tensor_cpp = new Tensor(6,user_coordinates,"Euler_ZXZ");
  const auto general_coordinates_zxz = t_tensor_cpp->UserCoordinatesToGeneralCoordinates();
  const std::vector<double> ar_cpp = Tensor::Ar(general_coordinates_zxz);

  for (int i=0;i<3;++i) {
    std::cout<<"ar: "<<i<<" -- "<<ar_test.at(i)<<" ; "<<ar_cpp.at(i)<<std::endl;
    EXPECT_EQ(true, std::Norm(ar_cpp,ar_test)<pow(10,-8));
  }



  infile.close();
}

TEST(Tensor, T_Tensor) {

  std::vector<double> from_text(9646,0);
  std::ifstream infile("../../../testing_datafiles/t_tensors/t_tensors_test.txt",std::ios_base::in);


  if (!infile) {
    std::cerr << " * Failed to open file "  << std::endl;
    if (infile.fail()) {
      // Print a more detailed error message using
      // strerror
      std::cerr << "Error details: " << strerror(errno)
           << std::endl;
    }
  }

  AuxiliarFunctions<double>::readArray(infile, from_text.data(), 9646);

  const auto coordinates_zxz =  std::vector<double>(from_text.begin(),from_text.begin()+6);
  const auto t_tensor_file =  std::vector<double>(from_text.begin()+6,from_text.end());

  auto *t_tensor_cpp = new Tensor(6, coordinates, "Euler_ZXZ");
  const std::vector<double> t_cpp = t_tensor_cpp->CalculateTensor(15);

  std::cout<<"T_Test: "<<t_tensor_file.at(0)<<" ; "<<t_cpp.at(0)<<std::endl;
  for (int i=0;i<1;i++) {
    EXPECT_EQ(t_tensor_file.at(i), t_cpp.at(i));
  }


  infile.close();
}


int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
