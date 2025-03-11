#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include <ALFI/poly.h>
#include <ALFI/dist.h>

#define POLYNOMIAL alfi::poly::lagrange
#define POLYNOMIAL_VALS alfi::poly::lagrange_vals
#include "test_case.h"

TEST(LagrangeTest, Basic) {
	test_case({0, 1, 2}, {0, 0, 0});
	test_case({0, 1, 2}, {1, 1, 1});
	test_case({0, 1, 2}, {1, 2, 3});
	test_case({0, 1, 2, 3, 5}, {3, 1, 2, 5, -1});
}

TEST(LagrangeTest, Functions) {
	test_case(
		alfi::dist::uniform(3, 0.0, 20.0),
		alfi::dist::uniform(10, 0.0, 20.0),
		[](double x) { return x * x; });
	test_case(
		alfi::dist::chebyshev(21, -2.0, 2.0),
		alfi::dist::uniform(10, -2.0, 2.0),
		[](double x) { return std::exp(x); },
		1e-6);
	test_case(
		alfi::dist::chebyshev_2(31, -2*M_PI, 2*M_PI),
		alfi::dist::uniform(10, -2*M_PI, 2*M_PI),
		[](double x) { return std::sin(x) + std::cos(x); },
		2e-2);
}