
//
// Created by albpl on 3/19/2025.
//

#include <gtest/gtest.h>


int Factorials(int n) {
  int result = 1;
  for (int i = 1; i <= n; i++) {
    result *= i;

  }
  return result;
}


// Demonstrate some basic assertions.
TEST(FactorialTest, HandlesZeroInput) {
  EXPECT_EQ(Factorials(0), 1);
}

// Tests factorial of positive numbers.
TEST(FactorialTest, HandlesPositiveInput) {
  EXPECT_EQ(Factorials(1), 1);
  EXPECT_EQ(Factorials(2), 2);
  EXPECT_EQ(Factorials(3), 6);
  EXPECT_EQ(Factorials(8), 40320);
}



int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
