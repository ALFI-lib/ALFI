#include <gtest/gtest.h>

#include <ALFI.h>

TEST(StepSplineTest, General) {
	const std::vector<double> X = {0, 1, 2, 3};
	const std::vector<double> Y = {5, 2, 10, 4};
	const auto left_spline = alfi::spline::StepSpline(X, Y, alfi::spline::StepSpline<>::Type::LEFT);
	const auto right_spline = alfi::spline::StepSpline(X, Y, alfi::spline::StepSpline<>::Type::RIGHT);
	const auto middle_spline = alfi::spline::StepSpline(X, Y, alfi::spline::StepSpline<>::Type::MIDDLE);
	const auto epsilon = 0.0001;
	for (size_t i = 0; i < X.size() - 1; ++i) {
		// LEFT exact points
		EXPECT_EQ(left_spline(X[i]), Y[i]);
		EXPECT_EQ(left_spline(X[i+1]), Y[i+1]);
		// LEFT between points
		EXPECT_EQ(left_spline(X[i] + epsilon), Y[i]);
		EXPECT_EQ(left_spline((X[i]+X[i+1])/2), Y[i]);
		EXPECT_EQ(left_spline(X[i+1] - epsilon), Y[i]);

		// RIGHT exact points
		EXPECT_EQ(right_spline(X[i]), Y[i]);
		EXPECT_EQ(right_spline(X[i+1]), Y[i+1]);
		// RIGHT between points
		EXPECT_EQ(right_spline(X[i] + epsilon), Y[i+1]);
		EXPECT_EQ(right_spline((X[i]+X[i+1])/2), Y[i+1]);
		EXPECT_EQ(right_spline(X[i+1] - epsilon), Y[i+1]);

		// MIDDLE exact points
		EXPECT_EQ(middle_spline(X[i]), Y[i]);
		EXPECT_EQ(middle_spline(X[i+1]), Y[i+1]);
		// MIDDLE between points
		EXPECT_EQ(middle_spline(X[i] + epsilon), (Y[i]+Y[i+1])/2);
		EXPECT_EQ(middle_spline((X[i]+X[i+1])/2), (Y[i]+Y[i+1])/2);
		EXPECT_EQ(middle_spline(X[i+1] - epsilon), (Y[i]+Y[i+1])/2);
	}
	EXPECT_EQ(left_spline(X.front() - 1), Y.front());
	EXPECT_EQ(left_spline(X.back() + 1), Y.back());
	EXPECT_EQ(right_spline(X.front() - 1), Y.front());
	EXPECT_EQ(right_spline(X.back() + 1), Y.back());
	EXPECT_EQ(middle_spline(X.front() - 1), Y.front());
	EXPECT_EQ(middle_spline(X.back() + 1), Y.back());
}

TEST(StepSplineTest, Moving) {
	std::vector<double> X = {0, 1, 2, 3};
	std::vector<double> Y = {5, 2, 10, 4};
	auto spline = alfi::spline::StepSpline(std::move(X), std::move(Y));
	EXPECT_TRUE(X.empty());
	EXPECT_TRUE(Y.empty());
	EXPECT_FALSE(spline.X().empty());
	EXPECT_FALSE(spline.Y().empty());
	X = spline.X();
	Y = spline.Y();
	EXPECT_FALSE(X.empty());
	EXPECT_FALSE(Y.empty());
	EXPECT_FALSE(spline.X().empty());
	EXPECT_FALSE(spline.Y().empty());
	X = std::move(spline).X();
	Y = std::move(spline).Y();
	EXPECT_FALSE(X.empty());
	EXPECT_FALSE(Y.empty());
	EXPECT_TRUE(spline.X().empty());
	EXPECT_TRUE(spline.Y().empty());
}