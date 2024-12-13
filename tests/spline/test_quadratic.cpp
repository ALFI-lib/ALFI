#include <gtest/gtest.h>

#include <ALFI.h>

TEST(QuadraticSplineTest, General) {
	const std::vector<double> X = {0, 1, 2, 3};
	const std::vector<double> Y = {5, 2, 10, 4};
	const auto spline = alfi::spline::QuadraticSpline(X, Y);
	for (size_t i = 0; i < X.size() - 1; ++i) {
		EXPECT_EQ(spline(X[i]), Y[i]);
		EXPECT_EQ(spline(X[i+1]), Y[i+1]);
	}
}

TEST(QuadraticSplineTest, Moving) {
	std::vector<double> X = {0, 1, 2, 3};
	const std::vector<double> Y = {5, 2, 10, 4};
	auto spline = alfi::spline::QuadraticSpline(std::move(X), Y);
	EXPECT_TRUE(X.empty());
	EXPECT_FALSE(spline.X().empty());
	EXPECT_FALSE(spline.coeffs().empty());
	X = spline.X();
	auto coeffs = spline.coeffs();
	EXPECT_FALSE(X.empty());
	EXPECT_FALSE(coeffs.empty());
	EXPECT_FALSE(spline.X().empty());
	EXPECT_FALSE(spline.coeffs().empty());
	X = std::move(spline).X();
	coeffs = std::move(spline).coeffs();
	EXPECT_FALSE(X.empty());
	EXPECT_FALSE(coeffs.empty());
	EXPECT_TRUE(spline.X().empty());
	EXPECT_TRUE(spline.coeffs().empty());
}