
//
// Created by albpl on 3/19/2025.
//

#include <gtest/gtest.h>
#include "../Tensor.h"





const int system_dimension = 6;
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


int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
