/**
  run_test.cpp

  Purpose: Run unit tests.

  @author Thomas Caleb

  @version 2.0 05/09/2024
*/

#include "gtest/gtest.h"
#include "../include/soda.h"

using namespace std;


int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
