#include <gtest/gtest.h>

#include <ALFI.h>

TEST(DistributionsTest, Uniform) {

	static const auto test_case_uniform = [](size_t n, double a, double b, const std::vector<double>& expected) {
		const auto dist = alfi::dist::uniform<double>(n, a, b);
		for (size_t index = 0; index < expected.size(); index++) {
			EXPECT_DOUBLE_EQ(dist[index], expected[index]);
		}
	};

	test_case_uniform(0, 0, 1, {});
	test_case_uniform(1, 0, 1, {0.5});
	test_case_uniform(2, 0, 1, {0, 1});
	test_case_uniform(3, 0, 1, {0, 0.5, 1});
	test_case_uniform(4, 0, 1, {0, 0.33333333333333333, 0.66666666666666667, 1});
	test_case_uniform(5, 0, 1, {0, 0.25, 0.5, 0.75, 1});
	test_case_uniform(11, 0, 1, {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0});

	test_case_uniform(0, -10, 10, {});
	test_case_uniform(1, -10, 10, {0});
	test_case_uniform(2, -10, 10, {-10, 10});
	test_case_uniform(3, -10, 10, {-10, 0, 10});
	test_case_uniform(4, -10, 10, {-10, -3.3333333333333333, 3.3333333333333333, 10});
	test_case_uniform(5, -10, 10, {-10, -5, 0, 5, 10});
	test_case_uniform(6, -10, 10, {-10, -6, -2, 2, 6, 10});
}

TEST(DistributionsTest, Chebyshev) {

	static const auto test_case_chebyshev = [](size_t n, double a, double b, const std::vector<double>& expected) {
		const auto dist = alfi::dist::chebyshev<double>(n, a, b);
		for (size_t index = 0; index < expected.size(); index++) {
			EXPECT_DOUBLE_EQ(dist[index], expected[index]);
		}
	};

	test_case_chebyshev(0, 0, 1, {});
	test_case_chebyshev(1, 0, 1, {0.5});
	test_case_chebyshev(2, 0, 1, {0.14644660940672627, 0.85355339059327373});

	{{
		const auto chebyshev = alfi::dist::chebyshev<double>(3, 0, 1);
		const auto uniform = alfi::dist::uniform<double>(3, 0, 1);
		EXPECT_LT(uniform[0], chebyshev[0]);
		EXPECT_EQ(chebyshev[1], 0.5);
		EXPECT_GT(uniform[2], chebyshev[2]);
	}}
	{{
		const auto chebyshev = alfi::dist::chebyshev<double>(4, 0, 1);
		const auto uniform = alfi::dist::uniform<double>(4, 0, 1);
		EXPECT_LT(uniform[0], chebyshev[0]);
		EXPECT_GT(uniform[1], chebyshev[1]);
		EXPECT_LT(uniform[2], chebyshev[2]);
		EXPECT_GT(uniform[3], chebyshev[3]);
	}}
}

TEST(DistributionsTest, ChebyshevStretched) {

	static const auto test_case_chebyshev_stretched = [](size_t n, double a, double b, const std::vector<double>& expected) {
		const auto dist = alfi::dist::chebyshev_stretched<double>(n, a, b);
		for (size_t index = 0; index < expected.size(); index++) {
			EXPECT_DOUBLE_EQ(dist[index], expected[index]);
		}
	};

	test_case_chebyshev_stretched(0, 0, 1, {});
	test_case_chebyshev_stretched(1, 0, 1, {0.5});
	test_case_chebyshev_stretched(2, 0, 1, {0, 1});
	test_case_chebyshev_stretched(3, 0, 1, {0, 0.5, 1});
}

TEST(DistributionsTest, ChebyshevEllipse) {

	static const auto test_case_chebyshev_ellipse = [](size_t n, double a, double b, double ratio, const std::vector<double>& expected) {
		const auto dist = alfi::dist::chebyshev_ellipse<double>(n, a, b, ratio);
		for (size_t index = 0; index < expected.size(); index++) {
			EXPECT_DOUBLE_EQ(dist[index], expected[index]);
		}
	};

	test_case_chebyshev_ellipse(0, 0, 1, 1, {});
	test_case_chebyshev_ellipse(1, 0, 1, 1, {0.5});

	{{
		const size_t n = 5;
		const double a = 0, b = 1, ratio = 0.5;
		const auto chebyshev_ellipse = alfi::dist::chebyshev_ellipse<double>(n, a, b, ratio);
		const auto chebyshev = alfi::dist::chebyshev<double>(n, a, b);
		EXPECT_GT(chebyshev_ellipse[0], chebyshev[0]);
		EXPECT_GT(chebyshev_ellipse[1], chebyshev[1]);
		EXPECT_EQ(chebyshev_ellipse[2], chebyshev[2]);
		EXPECT_LT(chebyshev_ellipse[3], chebyshev[3]);
		EXPECT_LT(chebyshev_ellipse[4], chebyshev[4]);
	}}
	{{
		const size_t n = 5;
		const double a = 0, b = 1, ratio = 2;
		const auto chebyshev_ellipse = alfi::dist::chebyshev_ellipse<double>(n, a, b, ratio);
		const auto chebyshev = alfi::dist::chebyshev<double>(n, a, b);
		EXPECT_LT(chebyshev_ellipse[0], chebyshev[0]);
		EXPECT_LT(chebyshev_ellipse[1], chebyshev[1]);
		EXPECT_EQ(chebyshev_ellipse[2], chebyshev[2]);
		EXPECT_GT(chebyshev_ellipse[3], chebyshev[3]);
		EXPECT_GT(chebyshev_ellipse[4], chebyshev[4]);
	}}
}

TEST(DistributionsTest, ChebyshevEllipseStretched) {

	static const auto test_case_chebyshev_ellipse_stretched = [](size_t n, double a, double b, double ratio, const std::vector<double>& expected) {
		const auto dist = alfi::dist::chebyshev_ellipse_stretched<double>(n, a, b, ratio);
		for (size_t index = 0; index < expected.size(); index++) {
			EXPECT_DOUBLE_EQ(dist[index], expected[index]);
		}
	};

	test_case_chebyshev_ellipse_stretched(0, 0, 1, 1, {});
	test_case_chebyshev_ellipse_stretched(1, 0, 1, 1, {0.5});
	test_case_chebyshev_ellipse_stretched(2, 0, 1, 1, {0, 1});
	test_case_chebyshev_ellipse_stretched(3, 0, 1, 1, {0, 0.5, 1});
}

TEST(DistributionsTest, CircleProjection) {

	static const auto test_case_circle_projection = [](size_t n, double a, double b, const std::vector<double>& expected) {
		const auto dist = alfi::dist::circle_proj<double>(n, a, b);
		for (size_t index = 0; index < expected.size(); index++) {
			EXPECT_DOUBLE_EQ(dist[index], expected[index]);
		}
	};

	test_case_circle_projection(0, 0, 1, {});
	test_case_circle_projection(1, 0, 1, {0.5});
	test_case_circle_projection(2, 0, 1, {0, 1});
	test_case_circle_projection(3, 0, 1, {0, 0.5, 1});
}

TEST(DistributionsTest, EllipseProjection) {

	static const auto test_case_ellipse_projection = [](size_t n, double a, double b, double B, const std::vector<double>& expected) {
		const auto dist = alfi::dist::ellipse_proj<double>(n, a, b, B);
		for (size_t index = 0; index < expected.size(); index++) {
			EXPECT_DOUBLE_EQ(dist[index], expected[index]);
		}
	};

	test_case_ellipse_projection(0, 0, 1, 1, {});
	test_case_ellipse_projection(1, 0, 1, 1, {0.5});
	test_case_ellipse_projection(2, 0, 1, 1, {0, 1});
	test_case_ellipse_projection(3, 0, 1, 1, {0, 0.5, 1});
}

TEST(DistributionsTest, Sigmoid) {

	static const auto test_case_sigmoid = [](size_t n, double a, double b, double steepness, const std::vector<double>& expected) {
		const auto dist = alfi::dist::sigmoid<double>(n, a, b, steepness);
		for (size_t index = 0; index < expected.size(); index++) {
			EXPECT_DOUBLE_EQ(dist[index], expected[index]);
		}
	};

	test_case_sigmoid(0, 0, 1, 1, {});
	test_case_sigmoid(1, 0, 1, 1, {0.5});
}

TEST(DistributionsTest, SigmoidStretched) {

	static const auto test_case_sigmoid_stretched = [](size_t n, double a, double b, double steepness, const std::vector<double>& expected) {
		const auto dist = alfi::dist::sigmoid_stretched<double>(n, a, b, steepness);
		for (size_t index = 0; index < expected.size(); index++) {
			EXPECT_DOUBLE_EQ(dist[index], expected[index]);
		}
	};

	test_case_sigmoid_stretched(0, 0, 1, 1, {});
	test_case_sigmoid_stretched(1, 0, 1, 1, {0.5});
	test_case_sigmoid_stretched(2, 0, 1, 1, {0, 1});
	test_case_sigmoid_stretched(3, 0, 1, 1, {0, 0.5, 1});
}

TEST(DistributionsTest, Erf) {

	static const auto test_case_erf = [](size_t n, double a, double b, double steepness, const std::vector<double>& expected) {
		const auto dist = alfi::dist::erf<double>(n, a, b, steepness);
		for (size_t index = 0; index < expected.size(); index++) {
			EXPECT_DOUBLE_EQ(dist[index], expected[index]);
		}
	};

	test_case_erf(0, 0, 1, 1, {});
	test_case_erf(1, 0, 1, 1, {0.5});
}

TEST(DistributionsTest, ErfStretched) {

	static const auto test_case_erf_stretched = [](size_t n, double a, double b, double steepness, const std::vector<double>& expected) {
		const auto dist = alfi::dist::erf_stretched<double>(n, a, b, steepness);
		for (size_t index = 0; index < expected.size(); index++) {
			EXPECT_DOUBLE_EQ(dist[index], expected[index]);
		}
	};

	test_case_erf_stretched(0, 0, 1, 1, {});
	test_case_erf_stretched(1, 0, 1, 1, {0.5});
	test_case_erf_stretched(2, 0, 1, 1, {0, 1});
	test_case_erf_stretched(3, 0, 1, 1, {0, 0.5, 1});
}