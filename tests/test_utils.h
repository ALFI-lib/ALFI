#pragma once

#include <gtest/gtest.h>

inline void expect_eq(const std::vector<double>& arr1, const std::vector<double>& arr2, double epsilon = 1e-12) {
	ASSERT_EQ(arr1.size(), arr2.size());
	for (std::size_t i = 0; i < arr1.size(); ++i) {
		EXPECT_NEAR(arr1[i], arr2[i], epsilon) << "i = " << i;
	}
}