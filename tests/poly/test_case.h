#pragma once

#include <vector>

#include "../test_utils.h"

inline void test_case(const std::vector<double>& X, const std::vector<double>& Y,
					  const std::vector<double>& expected_coeffs,
					  const std::vector<double>& xx, const std::vector<double>& yy,
					  double epsilon = 1e-9) {
	const auto P = POLYNOMIAL(X, Y);
	expect_eq(P, expected_coeffs);
	const auto values = alfi::poly::val(P, xx);
	expect_eq(values, yy, epsilon);
	const auto vals = POLYNOMIAL_VALS(X, Y, xx);
	expect_eq(vals, yy, epsilon);
}

inline void test_case(const std::vector<double>& X, const std::vector<double>& Y,
					  const std::vector<double>& xx, const std::vector<double>& yy,
					  double epsilon = 1e-9) {
	const auto P = POLYNOMIAL(X, Y);
	const auto values = alfi::poly::val(P, xx);
	expect_eq(values, yy, epsilon);
	const auto vals = POLYNOMIAL_VALS(X, Y, xx);
	expect_eq(vals, yy, epsilon);
}

inline void test_case(const std::vector<double>& X, const std::vector<double>& Y,
					  const std::vector<double>& expected_coeffs, double epsilon = 1e-9) {
	test_case(X, Y, expected_coeffs, X, Y, epsilon);
}

inline void test_case(const std::vector<double>& X, const std::vector<double>& Y, double epsilon = 1e-9) {
	test_case(X, Y, X, Y, epsilon);
}

inline void test_case(const std::vector<double>& X, const std::vector<double>& xx,
					  const std::function<double(double)>& f, double epsilon = 1e-9) {
	std::vector<double> Y(X.size());
	std::transform(X.begin(), X.end(), Y.begin(), f);
	std::vector<double> yy(xx.size());
	std::transform(xx.begin(), xx.end(), yy.begin(), f);
	test_case(X, Y, xx, yy, epsilon);
}