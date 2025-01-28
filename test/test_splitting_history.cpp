/**
	test_splitting_history.cpp

	Purpose: Test of the implementation of the SplittingHistory class.

	@author Thomas Caleb

	@version 1.0 23/01/2025
*/

#include "gtest/gtest.h"
#include "../include/soda.h"

using namespace DACE;
using namespace std;

TEST(TestSplittingHistory, ConstructorGetSet) {
    // Default constructor
    SplittingHistory history_default;
    EXPECT_EQ(history_default.size(), 0);

    // Copy constructor
    SplittingHistory history_to_copy;
    history_to_copy.push_back({1, 100});
    history_to_copy.push_back({2, -200});
    SplittingHistory history_copy(history_to_copy);
    EXPECT_EQ(history_copy.size(), history_to_copy.size());
    for (size_t i = 0; i < history_copy.size(); ++i) {
        EXPECT_EQ(history_copy[i], history_to_copy[i]);
    }

    // Copy assignment operator
    SplittingHistory history_assigned;
    history_assigned = history_to_copy;
    EXPECT_EQ(history_assigned.size(), history_to_copy.size());
    for (size_t i = 0; i < history_assigned.size(); ++i) {
        EXPECT_EQ(history_assigned[i], history_to_copy[i]);
    }
}
TEST(TestSplittingHistory, Alpha) {
    // Empty
    SplittingHistory history;
    EXPECT_EQ(history.alpha(), 1.0);

    // Add values
    history.push_back({1, 1});
    EXPECT_EQ(history.alpha(), ALPHA_1_GMM);


    history.push_back({2, -1});
    history.push_back({2, 0});
    EXPECT_NEAR(history.alpha(), ALPHA_0_GMM*ALPHA_1_GMM*ALPHA_1_GMM, EPS);
}
TEST(TestSplittingHistory, CanMerge) {
    SplittingHistory sh1, sh2, sh3;

    // Populate sh1 and sh2 with identical elements
    sh1.push_back({1, 2});
    sh1.push_back({3, 4});
    sh2.push_back({1, 2});
    sh2.push_back({3, 5});  // Different second element in the last pair

    // Populate sh3 with different size
    sh3.push_back({1, 2});

    // Test cases
    EXPECT_FALSE(sh1.can_merge(sh1));  // Should not be able to merge with itself
    EXPECT_TRUE(sh1.can_merge(sh2)); // Last element's second value is different
    EXPECT_FALSE(sh1.can_merge(sh3)); // Different sizes
}
TEST(TestSplittingHistory, StreamOperators) {
    // Test the << operator
    SplittingHistory history;
    history.push_back({1, 100});
    history.push_back({2, -200});

    std::ostringstream oss;
    oss << history;
    std::string expected_output = "1, 100; 2, -200; ";
    EXPECT_EQ(oss.str(), expected_output);

    // Test the >> operator
    std::istringstream iss("3, 300; 4, -400; ");
    SplittingHistory history_input;
    iss >> history_input;

    EXPECT_EQ(history_input.size(), 2);
    EXPECT_EQ(history_input[0], std::make_pair(3u, 300));
    EXPECT_EQ(history_input[1], std::make_pair(4u, -400));
}