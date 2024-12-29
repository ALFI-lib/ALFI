#include <gtest/gtest.h>

#include <ALFI/spline/linear.h>

TEST(LinearSplineTest, General) {
	const std::vector<double> X = {0, 1, 2, 3};
	const std::vector<double> Y = {5, 2, 10, 4};
	const auto spline = alfi::spline::LinearSpline(X, Y);
	const auto epsilon = 0.0001;
	for (size_t i = 0; i < X.size() - 1; ++i) {
		// exact points
		EXPECT_EQ(spline(X[i]), Y[i]);
		EXPECT_EQ(spline(X[i+1]), Y[i+1]);
		EXPECT_EQ(spline(X[i]), Y[i]);
		EXPECT_EQ(spline(X[i+1]), Y[i+1]);
		// between points
		EXPECT_DOUBLE_EQ(spline(X[i] + epsilon), Y[i] + epsilon * (Y[i+1] - Y[i]));
		EXPECT_DOUBLE_EQ(spline((X[i]+X[i+1])/2), (Y[i] + Y[i+1])/2);
		EXPECT_DOUBLE_EQ(spline(X[i+1] - epsilon), Y[i+1] - epsilon * (Y[i+1] - Y[i]));
	}
	EXPECT_EQ(spline(X[0] - 1), Y[0] - 1 * (Y[1] - Y[0]));
	EXPECT_EQ(spline(X[X.size()-1] + 1), Y[Y.size()-1] + 1 * (Y[Y.size()-1] - Y[Y.size()-2]));
}

TEST(LinearSplineTest, Moving) {
	std::vector<double> X = {0, 1, 2, 3};
	const std::vector<double> Y = {5, 2, 10, 4};
	auto spline = alfi::spline::LinearSpline(std::move(X), Y);
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