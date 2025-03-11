#include <gtest/gtest.h>

#include <ALFI/dist.h>
#include <ALFI/spline/cubic.h>

TEST(CubicSplineTest, General) {
	const std::vector<double> X = {0, 1, 2, 3};
	const std::vector<double> Y = {5, 2, 10, 4};
	const auto spline = alfi::spline::CubicSpline(X, Y);
	for (size_t i = 0; i < X.size() - 1; ++i) {
		EXPECT_EQ(spline(X[i]), Y[i]);
		EXPECT_EQ(spline(X[i+1]), Y[i+1]);
	}
}

TEST(CubicSplineTest, BasicConditions) {
	const size_t n = 10;
	const double a = -1, b = 1;
	const auto X = alfi::dist::uniform(n, a, b);
	const auto dX = alfi::util::arrays::diff(X);
	const auto Y = alfi::dist::chebyshev_2(n, a, b);
	auto spline = alfi::spline::CubicSpline<>(X, Y, alfi::spline::CubicSpline<>::Types::Natural{});
	EXPECT_NEAR(0, 2*spline.coeffs()[1], 1e-17);
	EXPECT_NEAR(0, 2*spline.coeffs()[4*8+1] + 6*spline.coeffs()[4*8]*dX[8], 1e-15);
	spline.construct(X, Y, alfi::spline::CubicSpline<>::Types::NotAKnot{});
	EXPECT_NEAR(spline.coeffs()[0], spline.coeffs()[4], 1e-15);
	EXPECT_NEAR(spline.coeffs()[4*7], spline.coeffs()[4*8], 1e-14);
	spline.construct(X, Y, alfi::spline::CubicSpline<>::Types::Periodic{});
	EXPECT_NEAR(spline.coeffs()[2], spline.coeffs()[4*8+2] + 2*spline.coeffs()[4*8+1]*dX[8] + 3*spline.coeffs()[4*8]*dX[8]*dX[8], 1e-16);
	EXPECT_NEAR(2*spline.coeffs()[1], 2*spline.coeffs()[4*8+1] + 6*spline.coeffs()[4*8]*dX[8], 1e-15);
	spline.construct(X, Y, alfi::spline::CubicSpline<>::Types::ParabolicEnds{});
	EXPECT_NEAR(0, spline.coeffs()[0], 1e-17);
	EXPECT_NEAR(0, spline.coeffs()[4*8], 1e-17);
	spline.construct(X, Y, alfi::spline::CubicSpline<>::Types::Clamped{-10, 10});
	EXPECT_NEAR(-10, spline.coeffs()[2], 1e-14);
	EXPECT_NEAR(10, spline.coeffs()[4*8+2] + 2*spline.coeffs()[4*8+1]*dX[8] + 3*spline.coeffs()[4*8]*dX[8]*dX[8], 1e-14);
	spline.construct(X, Y, alfi::spline::CubicSpline<>::Types::FixedSecond{-10, 10});
	EXPECT_NEAR(-10, 2*spline.coeffs()[1], 1e-17);
	EXPECT_NEAR(10, 2*spline.coeffs()[4*8+1] + 6*spline.coeffs()[4*8]*dX[8], 1e-17);
	spline.construct(X, Y, alfi::spline::CubicSpline<>::Types::FixedThird{-10, 10});
	EXPECT_NEAR(-10, 6*spline.coeffs()[0], 1e-17);
	EXPECT_NEAR(10, 6*spline.coeffs()[4*8], 1e-17);
	spline.construct(X, Y, alfi::spline::CubicSpline<>::Types::NotAKnotStart{});
	EXPECT_NEAR(spline.coeffs()[0], spline.coeffs()[4], 1e-17);
	EXPECT_NEAR(spline.coeffs()[4], spline.coeffs()[8], 1e-17);
	spline.construct(X, Y, alfi::spline::CubicSpline<>::Types::NotAKnotEnd{});
	EXPECT_NEAR(spline.coeffs()[4*6], spline.coeffs()[4*7], 1e-14);
	EXPECT_NEAR(spline.coeffs()[4*7], spline.coeffs()[4*8], 1e-15);
}

TEST(CubicSplineTest, CustomConditions) {
	const size_t n = 10;
	const double a = -1, b = 1;
	const auto X = alfi::dist::uniform(n, a, b);
	const auto dX = alfi::util::arrays::diff(X);
	const auto Y = alfi::dist::chebyshev_2(n, a, b);
	auto spline = alfi::spline::CubicSpline<>();
	// one point
	for (size_t i = 0; i < n; ++i) {
		spline.construct(X, Y, alfi::spline::CubicSpline<>::Types::Custom{
			alfi::spline::CubicSpline<>::Conditions::Clamped{i, 10},
			alfi::spline::CubicSpline<>::Conditions::FixedSecond{i, 10}});
		if (i < n - 1) {
			EXPECT_NEAR(10, spline.coeffs()[4*i+2], 1e-17);
			EXPECT_NEAR(10, 2*spline.coeffs()[4*i+1], 1e-17);
		}
		if (i > 0) {
			EXPECT_NEAR(10, spline.coeffs()[4*(i-1)+2] + 2*spline.coeffs()[4*(i-1)+1]*dX[(i-1)] + 3*spline.coeffs()[4*(i-1)]*dX[(i-1)]*dX[(i-1)], 1e-13);
			EXPECT_NEAR(10, 2*spline.coeffs()[4*(i-1)+1] + 6*spline.coeffs()[4*(i-1)]*dX[(i-1)], 1e-13);
		}
	}
	// two adjacent points
	for (size_t i = 0, j = 1; j < n; ++i, ++j) {
		spline.construct(X, Y, alfi::spline::CubicSpline<>::Types::Custom{
			alfi::spline::CubicSpline<>::Conditions::Clamped{i, 10},
			alfi::spline::CubicSpline<>::Conditions::FixedSecond{j, 10}});
		EXPECT_NEAR(10, spline.coeffs()[4*i+2], 1e-14);
		if (i > 0) {
			EXPECT_NEAR(10, spline.coeffs()[4*(i-1)+2] + 2*spline.coeffs()[4*(i-1)+1]*dX[(i-1)] + 3*spline.coeffs()[4*(i-1)]*dX[(i-1)]*dX[(i-1)], 1e-13);
		}
		if (j < n - 1) {
			EXPECT_NEAR(10, 2*spline.coeffs()[4*j+1], 1e-13);
		}
		EXPECT_NEAR(10, 2*spline.coeffs()[4*(j-1)+1] + 6*spline.coeffs()[4*(j-1)]*dX[(j-1)], 1e-13);
		spline.construct(X, Y, alfi::spline::CubicSpline<>::Types::Custom{
			alfi::spline::CubicSpline<>::Conditions::FixedThird{i, 10},
			alfi::spline::CubicSpline<>::Conditions::Clamped{j, 10}});
		EXPECT_NEAR(10, 6*spline.coeffs()[4*i], 1e-17);
		if (j < n - 1) {
			EXPECT_NEAR(10, spline.coeffs()[4*j+2], 1e-14);
		}
		EXPECT_NEAR(10, spline.coeffs()[4*(j-1)+2] + 2*spline.coeffs()[4*(j-1)+1]*dX[(j-1)] + 3*spline.coeffs()[4*(j-1)]*dX[(j-1)]*dX[(j-1)], 1e-13);
	}
	// various cases
	spline.construct(X, Y, alfi::spline::CubicSpline<>::Types::Custom{
		alfi::spline::CubicSpline<>::Conditions::Clamped{2, 10},
		alfi::spline::CubicSpline<>::Conditions::FixedSecond{5, 10}});
	EXPECT_NEAR(10, spline.coeffs()[4*2+2], 1e-14);
	EXPECT_NEAR(10, spline.coeffs()[4*1+2] + 2*spline.coeffs()[4*1+1]*dX[1] + 3*spline.coeffs()[4*1]*dX[1]*dX[1], 1e-13);
	EXPECT_NEAR(10, 2*spline.coeffs()[4*5+1], 1e-17);
	EXPECT_NEAR(10, 2*spline.coeffs()[4*4+1] + 6*spline.coeffs()[4*4]*dX[4], 1e-17);
	spline.construct(X, Y, alfi::spline::CubicSpline<>::Types::Custom{
		alfi::spline::CubicSpline<>::Conditions::NotAKnot{2},
		alfi::spline::CubicSpline<>::Conditions::FixedSecond{5, 10}});
	EXPECT_NEAR(spline.coeffs()[4*1], spline.coeffs()[4*2], 1e-17);
	EXPECT_NEAR(10, 2*spline.coeffs()[4*5+1], 1e-14);
	EXPECT_NEAR(10, 2*spline.coeffs()[4*4+1] + 6*spline.coeffs()[4*4]*dX[4], 1e-14);
	spline.construct(X, Y, alfi::spline::CubicSpline<>::Types::Custom{
		alfi::spline::CubicSpline<>::Conditions::NotAKnot{5},
		alfi::spline::CubicSpline<>::Conditions::FixedThird{5, 10}});
	EXPECT_NEAR(spline.coeffs()[4*4], spline.coeffs()[4*5], 1e-14);
	EXPECT_NEAR(10, 6*spline.coeffs()[4*4], 1e-14);
	EXPECT_NEAR(10, 6*spline.coeffs()[4*5], 1e-17);
	spline.construct(X, Y, alfi::spline::CubicSpline<>::Types::Custom{
		alfi::spline::CubicSpline<>::Conditions::FixedThird{2, -10},
		alfi::spline::CubicSpline<>::Conditions::FixedThird{5, 10}});
	EXPECT_NEAR(-10, 6*spline.coeffs()[4*2], 1e-17);
	EXPECT_NEAR(10, 6*spline.coeffs()[4*5], 1e-17);
	spline.construct(X, Y, alfi::spline::CubicSpline<>::Types::Custom{
		alfi::spline::CubicSpline<>::Conditions::FixedSecond{2, -10},
		alfi::spline::CubicSpline<>::Conditions::FixedThird{5, 10}});
	EXPECT_NEAR(-10, 2*spline.coeffs()[4*2+1], 1e-17);
	EXPECT_NEAR(-10, 2*spline.coeffs()[4*1+1] + 6*spline.coeffs()[4*1]*dX[1], 1e-17);
	EXPECT_NEAR(10, 6*spline.coeffs()[4*5], 1e-17);
	// same point (expect arithmetic mean)
	spline.construct(X, Y, alfi::spline::CubicSpline<>::Types::Custom{
		alfi::spline::CubicSpline<>::Conditions::Clamped{2, -10},
		alfi::spline::CubicSpline<>::Conditions::Clamped{2, 20}});
	EXPECT_NEAR(5, spline.coeffs()[4*2+2], 1e-17);
	EXPECT_NEAR(5, spline.coeffs()[4*1+2] + 2*spline.coeffs()[4*1+1]*dX[1] + 3*spline.coeffs()[4*1]*dX[1]*dX[1], 1e-14);
	spline.construct(X, Y, alfi::spline::CubicSpline<>::Types::Custom{
		alfi::spline::CubicSpline<>::Conditions::FixedSecond{2, -10},
		alfi::spline::CubicSpline<>::Conditions::FixedSecond{2, 20}});
	EXPECT_NEAR(5, 2*spline.coeffs()[4*2+1], 1e-17);
	EXPECT_NEAR(5, 2*spline.coeffs()[4*1+1] + 6*spline.coeffs()[4*1]*dX[1], 1e-17);
	spline.construct(X, Y, alfi::spline::CubicSpline<>::Types::Custom{
		alfi::spline::CubicSpline<>::Conditions::FixedThird{2, -10},
		alfi::spline::CubicSpline<>::Conditions::FixedThird{2, 20}});
	EXPECT_NEAR(5, 6*spline.coeffs()[4*2], 1e-17);
	spline.construct(X, Y, alfi::spline::CubicSpline<>::Types::Custom{
		alfi::spline::CubicSpline<>::Conditions::NotAKnot{2},
		alfi::spline::CubicSpline<>::Conditions::NotAKnot{2}});
	EXPECT_NEAR(spline.coeffs()[4*1], spline.coeffs()[4*2], 1e-15);
}

TEST(CubicSplineTest, Moving) {
	std::vector<double> X = {0, 1, 2, 3};
	const std::vector<double> Y = {5, 2, 10, 4};
	auto spline = alfi::spline::CubicSpline(std::move(X), Y);
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