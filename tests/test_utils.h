#pragma once

#include <gtest/gtest.h>

#include <toml++/toml.hpp>

#include <ALFI/util/numeric.h>

inline void expect_eq(const std::vector<double>& arr1, const std::vector<double>& arr2, double epsilon = 1e-12) {
	ASSERT_EQ(arr1.size(), arr2.size());
	for (std::size_t i = 0; i < arr1.size(); ++i) {
		EXPECT_TRUE(alfi::util::numeric::are_equal(arr1[i], arr2[i], epsilon))
			<< "i = " << i << "\narr1[i] = " << arr1[i] << "\narr2[i] = " << arr2[i] << "\nepsilon = " << epsilon;
	}
}

template <typename T>
std::vector<T> to_vector(const toml::array& array) {
	std::vector<T> result;
	result.reserve(array.size());
	for (const auto& value : array) {
		result.push_back(value.value<T>().value());
	}
	return result;
}