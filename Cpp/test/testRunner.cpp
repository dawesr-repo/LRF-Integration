
//
// Created by albpl on 3/19/2025.
//

#include <gtest/gtest.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <tuple>
#include <chrono>
#include "../AuxiliarFunctions.h"
#include "../Tensor.h"
#include "../PotentialEnergySurface.h"

constexpr int system_dimension = 6;
const std::vector <double> coordinates {10.00,10.0,20.0,30.0,40.0,50.0};
const std::string coordinate_format = "Euler_ZYZ";

auto* t_tensor = new Tensor(system_dimension,coordinates,coordinate_format);

TEST(Primitives_for_Tensor, Factorial) {

  EXPECT_EQ(0, t_tensor->Factorial(-1));
  EXPECT_EQ(1, t_tensor->Factorial(0));
  EXPECT_EQ(1, t_tensor->Factorial(1));
  EXPECT_EQ(24, t_tensor->Factorial(4));

}

TEST(Primitives_for_Tensor, FactorialNN) {

  EXPECT_EQ(0, t_tensor->FactorialNN(-1,1,1,1));
  EXPECT_EQ(0, t_tensor->FactorialNN(1,-1,1,1));
  EXPECT_EQ(0, t_tensor->FactorialNN(1,1,-1,1));
  EXPECT_EQ(0, t_tensor->FactorialNN(1,1,1,-1));
  EXPECT_EQ(0, t_tensor->FactorialNN(1,2,1,1));
  EXPECT_EQ(0, t_tensor->FactorialNN(1,1,1,2));
  EXPECT_EQ(1, t_tensor->FactorialNN(0,0,0,0));
  EXPECT_EQ(2, t_tensor->FactorialNN(1,1,1,1));

}

TEST(Primitives_for_Tensor, CoeffM) {

  EXPECT_EQ(1, t_tensor->CoeffM("z",0,"0"));
  EXPECT_EQ(1, t_tensor->CoeffM("z",1,"0"));
  EXPECT_EQ(1, t_tensor->CoeffM("z",2,"c"));
  EXPECT_EQ(1, t_tensor->CoeffM("z",3,"s"));

  EXPECT_EQ(0, t_tensor->CoeffM("x",0,"0"));
  EXPECT_EQ(sqrt(2.0), t_tensor->CoeffM("x",1,"c"));
  EXPECT_EQ(0, t_tensor->CoeffM("x",1,"s"));
  EXPECT_EQ(2, t_tensor->CoeffM("x",2,"c"));
  EXPECT_EQ(3, t_tensor->CoeffM("x",3,"s"));

  EXPECT_EQ(0, t_tensor->CoeffM("y",0,"0"));
  EXPECT_EQ(sqrt(2.0), t_tensor->CoeffM("y",1,"s"));
  EXPECT_EQ(0, t_tensor->CoeffM("y",1,"c"));
  EXPECT_EQ(-2, t_tensor->CoeffM("y",2,"c"));
  EXPECT_EQ(3, t_tensor->CoeffM("y",3,"s"));
}

TEST(Primitives_for_Tensor, GetComponent) {
  EXPECT_EQ(0, AuxiliarFunctions<int>::getComponents_v2(0,0,0,0));
  EXPECT_EQ(1, AuxiliarFunctions<int>::getComponents_v2(0,1,0,0));
  EXPECT_EQ(2, AuxiliarFunctions<int>::getComponents_v2(0,1,0,1));
  EXPECT_EQ(3, AuxiliarFunctions<int>::getComponents_v2(0,1,0,2));
  EXPECT_EQ(4, AuxiliarFunctions<int>::getComponents_v2(1,0,0,0));
  EXPECT_EQ(5, AuxiliarFunctions<int>::getComponents_v2(1,0,1,0));
  EXPECT_EQ(6, AuxiliarFunctions<int>::getComponents_v2(1,0,2,0));
}

TEST(Primitives_for_Tensor, NEta) {
  EXPECT_EQ(true, std::make_tuple(1,"0" ) == t_tensor->NEta("z",1,"0"));
  EXPECT_EQ(true, std::make_tuple(2,"c" ) == t_tensor->NEta("z",2,"c"));
  EXPECT_EQ(true, std::make_tuple(3,"s" ) == t_tensor->NEta("z",3,"s"));

  EXPECT_EQ(true, std::make_tuple(0,"0" ) == t_tensor->NEta("x",0,"0"));
  EXPECT_EQ(true, std::make_tuple(0,"0" ) == t_tensor->NEta("x",1,"c"));
  EXPECT_EQ(true, std::make_tuple(0,"0" ) == t_tensor->NEta("x",1,"s"));
  EXPECT_EQ(true, std::make_tuple(2,"c" ) == t_tensor->NEta("x",3,"c"));
  EXPECT_EQ(true, std::make_tuple(2,"s" ) == t_tensor->NEta("x",3,"s"));

  EXPECT_EQ(true, std::make_tuple(0,"0" ) == t_tensor->NEta("y",0,"0"));
  EXPECT_EQ(true, std::make_tuple(0,"0" ) == t_tensor->NEta("y",1,"c"));
  EXPECT_EQ(true, std::make_tuple(0,"0" ) == t_tensor->NEta("y",1,"s"));
  EXPECT_EQ(true, std::make_tuple(2,"s" ) == t_tensor->NEta("y",3,"c"));
  EXPECT_EQ(true, std::make_tuple(2,"c" ) == t_tensor->NEta("y",3,"s"));

}

TEST(Primitives_for_Tensor, GetSplittingComponent) {

  EXPECT_EQ("-1", t_tensor->GetSplittingComponent(-1));
  EXPECT_EQ("0", t_tensor->GetSplittingComponent(0));
  EXPECT_EQ("c", t_tensor->GetSplittingComponent(1));
  EXPECT_EQ("s", t_tensor->GetSplittingComponent(2));

}

TEST(Primitives_for_Tensor, GetTensorComponent) {

  EXPECT_EQ(-1, t_tensor->GetTensorComponent(-1,0,"0"));
  EXPECT_EQ(-1, t_tensor->GetTensorComponent(1,2,"0"));
  EXPECT_EQ(-1, t_tensor->GetTensorComponent(0,0,"c"));
  EXPECT_EQ(0, t_tensor->GetTensorComponent(4,0,"0"));
  EXPECT_EQ(1, t_tensor->GetTensorComponent(4,1,"c"));
  EXPECT_EQ(4, t_tensor->GetTensorComponent(4,2,"s"));

}

TEST(Primitives_for_Tensor, Ar_Components) {

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
  for (int i=0; i<100;++i) {
    AuxiliarFunctions<double>::readArray(infile, from_text.data(), 9);

    const auto user_coordinates =  std::vector<double>(from_text.begin(),from_text.begin()+6);
    auto ar_test = std::vector<double>(from_text.begin() + 6, from_text.end());


    const auto *t_tensor_cpp = new Tensor(6,user_coordinates,"Euler_ZXZ");
    const auto general_coordinates_zxz = t_tensor_cpp->UserCoordinatesToGeneralCoordinates();
    auto ar_cpp = Tensor::Ar(general_coordinates_zxz);

    EXPECT_EQ(true,AuxiliarFunctions<double>::areEqual(ar_cpp.data(),ar_test.data(), 3));
  }


  infile.close();
}

TEST(Primitives_for_Tensor,Br_Component) {
  std::vector<double> from_text(9,0);
  std::ifstream infile("../../../testing_datafiles/t_tensors/Br.txt",std::ios_base::in);


  if (!infile) {
    std::cerr << " * Failed to open file "  << std::endl;
    if (infile.fail()) {
      // Print a more detailed error message using
      // strerror
      std::cerr << "Error details: " << strerror(errno)
           << std::endl;
    }
  }
  for (int i=0; i<100;++i) {
    AuxiliarFunctions<double>::readArray(infile, from_text.data(), 9);

    const auto user_coordinates =  std::vector<double>(from_text.begin(),from_text.begin()+6);
    auto br_test = std::vector<double>(from_text.begin() + 6, from_text.end());


    const auto *t_tensor_cpp = new Tensor(6,user_coordinates,"Euler_ZXZ");
    const auto general_coordinates_zxz = t_tensor_cpp->UserCoordinatesToGeneralCoordinates();
    auto br_cpp = Tensor::Br(general_coordinates_zxz);

    EXPECT_EQ(true,AuxiliarFunctions<double>::areEqual(br_cpp.data(),br_test.data(), 3));
  }


  infile.close();
}

TEST(Primitives_for_Tensor,CC_Component) {
  std::vector<double> from_text(15,0);
  std::ifstream infile("../../../testing_datafiles/t_tensors/CC.txt");


  if (!infile) {
    std::cerr << " * Failed to open file "  << std::endl;
    if (infile.fail()) {
      // Print a more detailed error message using
      // strerror
      std::cerr << "Error details: " << strerror(errno)
           << std::endl;
    }
  }
  for (int i=0; i<100;++i) {
    AuxiliarFunctions<double>::readArray(infile, from_text.data(), 15);

    const auto user_coordinates =  std::vector<double>(from_text.begin(),from_text.begin()+6);
    auto cc_test = std::vector<double>(from_text.begin() + 6, from_text.end());


    const auto *t_tensor_cpp = new Tensor(6,user_coordinates,"Euler_ZXZ");
    const auto general_coordinates_zxz = t_tensor_cpp->UserCoordinatesToGeneralCoordinates();
    auto cc_cpp = Tensor::Cab(general_coordinates_zxz);

    EXPECT_EQ(true,AuxiliarFunctions<double>::areEqual(cc_cpp.data(),cc_test.data(), 4));
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
  for (int i=0; i<100;++i) {
    AuxiliarFunctions<double>::readArray(infile, from_text.data(), 9646);

    const auto coordinates_zxz =  std::vector<double>(from_text.begin(),from_text.begin()+6);
    auto t_tensor_file =  std::vector<double>(from_text.begin()+6,from_text.end());

    auto *t_tensor_cpp = new Tensor(6, coordinates_zxz, "Euler_ZXZ");
    std::vector<double> t_cpp = t_tensor_cpp->CalculateTensor(15);
    EXPECT_EQ(true,AuxiliarFunctions<double>::areEqual(t_cpp.data(),t_tensor_file.data(), 9640));
  }
}

TEST(Interaction_Energy, C1_C1_ions) {
  std::vector<double> from_text(10,0);
  std::ifstream infile(
      "../../../testing_datafiles/t_tensors/PES_Components.txt",
      std::ios_base::in);
  const auto * pes1 = new PotentialEnergySurface("../../../testing_datafiles/coefficients/C1(1)_C1(1)_Coeff.txt");

  if (!infile) {
    std::cerr << " * Failed to open file "  << std::endl;
    if (infile.fail()) {
      // Print a more detailed error message using
      // strerror
      std::cerr << "Error details: " << strerror(errno)
           << std::endl;
    }
  }



  for (int i=0; i<3;++i) {
    AuxiliarFunctions<double>::readArray(infile, from_text.data(), 10);

    const auto coordinates_zxz =  std::vector<double>(from_text.begin(),from_text.begin()+6);
    auto interactions_from_file =  std::vector<double>(from_text.begin()+6,from_text.end());

    auto *t_tensor = new Tensor(6, coordinates_zxz , "Euler_ZXZ");
    const std::vector<double> t = t_tensor->CalculateTensor(15);
    const double r = coordinates_zxz.at(0);

    const double multipole_interaction =  pes1->MultipoleInteraction(r, t);
    const double induction_interaction =  pes1->InductionInteraction(r, t);
    const double dispersion_interaction =  pes1->DispersionInteraction(r, t);
    const double total_interaction = multipole_interaction+induction_interaction+dispersion_interaction;

    EXPECT_EQ(true,std::abs(multipole_interaction-interactions_from_file.at(0))<DBL_EPSILON);
    EXPECT_EQ(true,std::abs(induction_interaction-interactions_from_file.at(1))<DBL_EPSILON);
    EXPECT_EQ(true,std::abs(dispersion_interaction-interactions_from_file.at(2))<DBL_EPSILON);
    EXPECT_EQ(true,std::abs(total_interaction-interactions_from_file.at(3))<DBL_EPSILON);
  }
}

TEST(Profile,EvaluateLRF) {
  const std::string path_file_pes1 =
  "../../../testing_datafiles/coefficients/C1(1)_C1(1)_Coeff.txt";
  auto* pes1 = new PotentialEnergySurface(path_file_pes1);
  constexpr int system_dimension = 6;
  const std::vector<double> coordinates{10.27, 30.0, 20.0, 120.0, 0.0, 0.0};
  const std::string coordinate_format = "Euler_ZYZ";
  const int n_test = 100;
  double accumative_time {0.0};
  for (int i=0;i<n_test;++i) {
    auto start = std::chrono::high_resolution_clock::now();
    pes1->EvaluateLRF(system_dimension, coordinates, coordinate_format);
    auto stop = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    accumative_time = accumative_time + duration.count();
  }

  std::cout << "Average Elapsed(ms)= " << accumative_time/n_test<< std::endl;

}


TEST(Components,versions) {
  int cpn{0};
  for (int order = 1; order <= 15; ++order) {
    for (int lap = 0; lap <= order - 1; ++lap) {
      int lbp = order - lap - 1;
      for (int kap = 0; kap <= 2 * lap; ++kap) {
        for (int kbp = 0; kbp <= 2 * lbp; ++kbp) {
          const int cpn_v2 = AuxiliarFunctions<int>::getComponents_v2(lap, lbp, kap, kbp);
         EXPECT_EQ(cpn,cpn_v2 );
          cpn++;

        }
      }
    }
  }
}


int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
