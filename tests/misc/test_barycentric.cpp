#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include <ALFI/misc.h>

template <typename Number>
void test_case(const std::vector<Number>& X, const std::vector<Number>& Y,
			   const std::vector<Number>& xx,
			   alfi::dist::Type dist_type,
			   const std::vector<Number>& yy) {
	EXPECT_EQ(alfi::misc::barycentric(X, Y, xx, dist_type), yy);
}

template <typename Number>
void test_case(const std::vector<Number>& X, const std::vector<Number>& Y, const std::vector<Number>& xx, const std::vector<Number>& yy) {
	EXPECT_EQ(alfi::misc::barycentric(X, Y, xx), yy);
}

template <typename Number>
void test_case(const std::vector<Number>& X, const std::vector<Number>& Y, alfi::dist::Type dist_type) {
	EXPECT_EQ(alfi::misc::barycentric(X, Y, X, dist_type), Y);
}

template <typename Number>
void test_case(const std::vector<Number>& X, const std::vector<Number>& Y) {
	EXPECT_EQ(alfi::misc::barycentric(X, Y, X), Y);
}

TEST(BarycentricTest, Basic) {
	test_case<double>({0, 1, 2}, {0, 0, 0});
	test_case<double>({0, 1, 2}, {1, 1, 1});
	test_case<double>({0, 1, 2}, {1, 2, 3});
}