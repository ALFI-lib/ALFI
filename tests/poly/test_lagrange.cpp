#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include <ALFI.h>

template <typename Number>
void test_case(const std::vector<Number>& X, const std::vector<Number>& Y, Number x0,
			   const std::vector<Number>& expected_coeffs,
			   const std::vector<Number>& xx, const std::vector<Number>& yy) {
	const auto L = alfi::poly::lagrange(X, Y, x0);
	EXPECT_EQ(L, expected_coeffs);
	const auto values = alfi::poly::val(L, xx, x0);
	EXPECT_EQ(values, yy);
}

template <typename Number>
void test_case(const std::vector<Number>& X, const std::vector<Number>& Y,
			   const std::vector<Number>& expected_coeffs,
			   const std::vector<Number>& xx, const std::vector<Number>& yy) {
	const auto L = alfi::poly::lagrange(X, Y);
	EXPECT_EQ(L, expected_coeffs);
	const auto values = alfi::poly::val(L, xx);
	EXPECT_EQ(values, yy);
}

template <typename Number>
void test_case(const std::vector<Number>& X, const std::vector<Number>& Y, Number x0,
			   const std::vector<Number>& xx, const std::vector<Number>& yy) {
	const auto L = alfi::poly::lagrange(X, Y, x0);
	const auto values = alfi::poly::val(L, xx, x0);
	EXPECT_EQ(values, yy);
}

template <typename Number>
void test_case(const std::vector<Number>& X, const std::vector<Number>& Y,
			   const std::vector<Number>& xx, const std::vector<Number>& yy) {
	const auto L = alfi::poly::lagrange(X, Y);
	const auto values = alfi::poly::val(L, xx);
	EXPECT_EQ(values, yy);
}

template <typename Number>
void test_case(const std::vector<Number>& X, const std::vector<Number>& Y, Number x0, std::vector<Number> expected_coeffs) {
	test_case(X, Y, x0, expected_coeffs, X, Y);
}

template <typename Number>
void test_case(const std::vector<Number>& X, const std::vector<Number>& Y, Number x0) {
	test_case(X, Y, x0, X, Y);
}

template <typename Number>
void test_case(const std::vector<Number>& X, const std::vector<Number>& Y, std::vector<Number> expected_coeffs) {
	test_case(X, Y, expected_coeffs, X, Y);
}

template <typename Number>
void test_case(const std::vector<Number>& X, const std::vector<Number>& Y) {
	test_case(X, Y, X, Y);
}

TEST(LagrangeTest, Basic) {
	test_case<double>({0, 1, 2}, {0, 0, 0});
	test_case<double>({0, 1, 2}, {1, 1, 1});
	test_case<double>({0, 1, 2}, {1, 2, 3});

	test_case<double>({0, 1, 2}, {0, 0, 0}, 5);
	test_case<double>({0, 1, 2}, {1, 1, 1}, 5);
	test_case<double>({0, 1, 2}, {1, 2, 3}, 5);
}