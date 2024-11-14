#include <ALFI.h>

#include <gtest/gtest.h>

TEST(DistributionsTest, Iterators) {

	static const auto test_case_iterators = [](size_t n, double a, double b) {
		const alfi::dist::Uniform<double> dist = alfi::dist::Uniform<double>(n, a, b);
		size_t i = 0;
		auto dist_iter = dist.begin();
		while (i < n) {
			EXPECT_DOUBLE_EQ(dist[i], *dist_iter);
			++i;
			++dist_iter;
		}
		EXPECT_EQ(dist_iter, dist.end());
	};

	test_case_iterators(0, 0, 1);
	test_case_iterators(3, 0, 1);
	test_case_iterators(10, 0, 1);

	test_case_iterators(0, 10, 100);
	test_case_iterators(3, 100, 10);
	test_case_iterators(10, 0.1234, 4321.9876);
}

TEST(DistributionsTest, UniformDistribution) {

	static const auto test_case_uniform = [](size_t n, double a, double b, const std::vector<double>& expected) {
		const alfi::dist::Uniform<double> dist = alfi::dist::Uniform<double>(n, a, b);
		EXPECT_EQ(dist.type(), alfi::dist::Type::UNIFORM);
		EXPECT_EQ(dist.n(), n);
		EXPECT_EQ(dist.a(), a);
		EXPECT_EQ(dist.b(), b);
		EXPECT_EQ(dist.n(), expected.size());
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

TEST(DistributionsTest, ChebyshevDistribution) {

	static const auto test_case_chebyshev = [](size_t n, double a, double b, const std::vector<double>& expected) {
		const alfi::dist::Chebyshev<double> dist = alfi::dist::Chebyshev<double>(n, a, b);
		EXPECT_EQ(dist.type(), alfi::dist::Type::CHEBYSHEV);
		EXPECT_EQ(dist.n(), n);
		EXPECT_EQ(dist.a(), a);
		EXPECT_EQ(dist.b(), b);
		EXPECT_EQ(dist.n(), expected.size());
		for (size_t index = 0; index < expected.size(); index++) {
			EXPECT_DOUBLE_EQ(dist[index], expected[index]);
		}
	};

	test_case_chebyshev(0, 0, 1, {});
	test_case_chebyshev(1, 0, 1, {0.5});
	test_case_chebyshev(2, 0, 1, {0.14644660940672627, 0.85355339059327373});

	{{
		const alfi::dist::Chebyshev<double> chebyshev = alfi::dist::Chebyshev<double>(3, 0, 1);
		const alfi::dist::Uniform<double> uniform = alfi::dist::Uniform<double>(3, 0, 1);
		EXPECT_EQ(chebyshev.n(), 3);
		EXPECT_LT(uniform[0], chebyshev[0]);
		EXPECT_EQ(chebyshev[1], 0.5);
		EXPECT_GT(uniform[2], chebyshev[2]);
	}}
	{{
		const alfi::dist::Chebyshev<double> chebyshev = alfi::dist::Chebyshev<double>(4, 0, 1);
		const alfi::dist::Uniform<double> uniform = alfi::dist::Uniform<double>(4, 0, 1);
		EXPECT_EQ(chebyshev.n(), 4);
		EXPECT_LT(uniform[0], chebyshev[0]);
		EXPECT_GT(uniform[1], chebyshev[1]);
		EXPECT_LT(uniform[2], chebyshev[2]);
		EXPECT_GT(uniform[3], chebyshev[3]);
	}}
}

TEST(DistributionsTest, ChebyshevStretchedDistribution) {

	static const auto test_case_chebyshev_stretched = [](size_t n, double a, double b, const std::vector<double>& expected) {
		const alfi::dist::ChebyshevStretched<double> dist = alfi::dist::ChebyshevStretched<double>(n, a, b);
		EXPECT_EQ(dist.type(), alfi::dist::Type::CHEBYSHEV_STRETCHED);
		EXPECT_EQ(dist.n(), n);
		EXPECT_EQ(dist.a(), a);
		EXPECT_EQ(dist.b(), b);
		EXPECT_LE(dist.stretched_a(), a);
		EXPECT_GE(dist.stretched_b(), b);
		EXPECT_EQ(dist.n(), expected.size());
		for (size_t index = 0; index < expected.size(); index++) {
			EXPECT_DOUBLE_EQ(dist[index], expected[index]);
		}
	};

	test_case_chebyshev_stretched(0, 0, 1, {});
	test_case_chebyshev_stretched(1, 0, 1, {0.5});
	test_case_chebyshev_stretched(2, 0, 1, {0, 1});
	test_case_chebyshev_stretched(3, 0, 1, {0, 0.5, 1});
}

TEST(DistributionsTest, ChebyshevEllipseDistribution) {

	static const auto test_case_chebyshev_ellipse = [](size_t n, double a, double b, double B, const std::vector<double>& expected) {
		const alfi::dist::ChebyshevEllipse<double> dist = alfi::dist::ChebyshevEllipse<double>(n, a, b, B);
		EXPECT_EQ(dist.type(), alfi::dist::Type::CHEBYSHEV_ELLIPSE);
		EXPECT_EQ(dist.n(), n);
		EXPECT_EQ(dist.a(), a);
		EXPECT_EQ(dist.b(), b);
		EXPECT_EQ(dist.ratio(), B);
		EXPECT_EQ(dist.n(), expected.size());
		for (size_t index = 0; index < expected.size(); index++) {
			EXPECT_DOUBLE_EQ(dist[index], expected[index]);
		}
	};

	test_case_chebyshev_ellipse(0, 0, 1, 1, {});
	test_case_chebyshev_ellipse(1, 0, 1, 1, {0.5});

	{{
		const size_t n = 5;
		const double a = 0, b = 1, B = 0.5;
		const alfi::dist::ChebyshevEllipse<double> chebyshev_ellipse = alfi::dist::ChebyshevEllipse<double>(n, a, b, B);
		const alfi::dist::Chebyshev<double> chebyshev = alfi::dist::Chebyshev<double>(n, a, b);
		EXPECT_GT(chebyshev_ellipse[0], chebyshev[0]);
		EXPECT_GT(chebyshev_ellipse[1], chebyshev[1]);
		EXPECT_EQ(chebyshev_ellipse[2], chebyshev[2]);
		EXPECT_LT(chebyshev_ellipse[3], chebyshev[3]);
		EXPECT_LT(chebyshev_ellipse[4], chebyshev[4]);
	}}
	{{
		const size_t n = 5;
		const double a = 0, b = 1, B = 2;
		const alfi::dist::ChebyshevEllipse<double> chebyshev_ellipse = alfi::dist::ChebyshevEllipse<double>(n, a, b, B);
		const alfi::dist::Chebyshev<double> chebyshev = alfi::dist::Chebyshev<double>(n, a, b);
		EXPECT_LT(chebyshev_ellipse[0], chebyshev[0]);
		EXPECT_LT(chebyshev_ellipse[1], chebyshev[1]);
		EXPECT_EQ(chebyshev_ellipse[2], chebyshev[2]);
		EXPECT_GT(chebyshev_ellipse[3], chebyshev[3]);
		EXPECT_GT(chebyshev_ellipse[4], chebyshev[4]);
	}}
}

TEST(DistributionsTest, ChebyshevEllipseStretchedDistribution) {

	static const auto test_case_chebyshev_ellipse_stretched = [](size_t n, double a, double b, double B, const std::vector<double>& expected) {
		const alfi::dist::ChebyshevEllipseStretched<double> dist = alfi::dist::ChebyshevEllipseStretched<double>(n, a, b, B);
		EXPECT_EQ(dist.type(), alfi::dist::Type::CHEBYSHEV_ELLIPSE_STRETCHED);
		EXPECT_EQ(dist.n(), n);
		EXPECT_EQ(dist.a(), a);
		EXPECT_EQ(dist.b(), b);
		EXPECT_EQ(dist.ratio(), B);
		EXPECT_LE(dist.stretched_a(), a);
		EXPECT_GE(dist.stretched_b(), b);
		EXPECT_EQ(dist.n(), expected.size());
		for (size_t index = 0; index < expected.size(); index++) {
			EXPECT_DOUBLE_EQ(dist[index], expected[index]);
		}
	};

	test_case_chebyshev_ellipse_stretched(0, 0, 1, 1, {});
	test_case_chebyshev_ellipse_stretched(1, 0, 1, 1, {0.5});
	test_case_chebyshev_ellipse_stretched(2, 0, 1, 1, {0, 1});
	test_case_chebyshev_ellipse_stretched(3, 0, 1, 1, {0, 0.5, 1});
}

TEST(DistributionsTest, CircleProjectionDistribution) {

	static const auto test_case_circle_projection = [](size_t n, double a, double b, const std::vector<double>& expected) {
		const alfi::dist::CircleProj<double> dist = alfi::dist::CircleProj<double>(n, a, b);
		EXPECT_EQ(dist.type(), alfi::dist::Type::CIRCLE_PROJECTION);
		EXPECT_EQ(dist.n(), n);
		EXPECT_EQ(dist.a(), a);
		EXPECT_EQ(dist.b(), b);
		EXPECT_EQ(dist.n(), expected.size());
		for (size_t index = 0; index < expected.size(); index++) {
			EXPECT_DOUBLE_EQ(dist[index], expected[index]);
		}
	};

	test_case_circle_projection(0, 0, 1, {});
	test_case_circle_projection(1, 0, 1, {0.5});
	test_case_circle_projection(2, 0, 1, {0, 1});
	test_case_circle_projection(3, 0, 1, {0, 0.5, 1});
}

TEST(DistributionsTest, EllipseProjectionDistribution) {

	static const auto test_case_ellipse_projection = [](size_t n, double a, double b, double B, const std::vector<double>& expected) {
		const alfi::dist::EllipseProj<double> dist = alfi::dist::EllipseProj<double>(n, a, b, B);
		EXPECT_EQ(dist.type(), alfi::dist::Type::ELLIPSE_PROJECTION);
		EXPECT_EQ(dist.n(), n);
		EXPECT_EQ(dist.a(), a);
		EXPECT_EQ(dist.b(), b);
		EXPECT_EQ(dist.ratio(), B);
		EXPECT_EQ(dist.n(), expected.size());
		for (size_t index = 0; index < expected.size(); index++) {
			EXPECT_DOUBLE_EQ(dist[index], expected[index]);
		}
	};

	test_case_ellipse_projection(0, 0, 1, 1, {});
	test_case_ellipse_projection(1, 0, 1, 1, {0.5});
	test_case_ellipse_projection(2, 0, 1, 1, {0, 1});
	test_case_ellipse_projection(3, 0, 1, 1, {0, 0.5, 1});
}

TEST(DistributionsTest, SigmoidDistribution) {

	static const auto test_case_sigmoid = [](size_t n, double a, double b, double steepness, const std::vector<double>& expected) {
		const alfi::dist::Sigmoid<double> dist = alfi::dist::Sigmoid<double>(n, a, b, steepness);
		EXPECT_EQ(dist.type(), alfi::dist::Type::SIGMOID);
		EXPECT_EQ(dist.n(), n);
		EXPECT_EQ(dist.a(), a);
		EXPECT_EQ(dist.b(), b);
		EXPECT_EQ(dist.steepness(), steepness);
		EXPECT_EQ(dist.n(), expected.size());
		for (size_t index = 0; index < expected.size(); index++) {
			EXPECT_DOUBLE_EQ(dist[index], expected[index]);
		}
	};

	test_case_sigmoid(0, 0, 1, 1, {});
	test_case_sigmoid(1, 0, 1, 1, {0.5});
}

TEST(DistributionsTest, SigmoidStretchedDistribution) {

	static const auto test_case_sigmoid_stretched = [](size_t n, double a, double b, double steepness, const std::vector<double>& expected) {
		const alfi::dist::SigmoidStretched<double> dist = alfi::dist::SigmoidStretched<double>(n, a, b, steepness);
		EXPECT_EQ(dist.type(), alfi::dist::Type::SIGMOID_STRETCHED);
		EXPECT_EQ(dist.n(), n);
		EXPECT_EQ(dist.a(), a);
		EXPECT_EQ(dist.b(), b);
		EXPECT_EQ(dist.steepness(), steepness);
		EXPECT_LE(dist.stretched_a(), a);
		EXPECT_GE(dist.stretched_b(), b);
		EXPECT_EQ(dist.n(), expected.size());
		for (size_t index = 0; index < expected.size(); index++) {
			EXPECT_DOUBLE_EQ(dist[index], expected[index]);
		}
	};

	test_case_sigmoid_stretched(0, 0, 1, 1, {});
	test_case_sigmoid_stretched(1, 0, 1, 1, {0.5});
	test_case_sigmoid_stretched(2, 0, 1, 1, {0, 1});
	test_case_sigmoid_stretched(3, 0, 1, 1, {0, 0.5, 1});
}

TEST(DistributionsTest, ErfDistribution) {

	static const auto test_case_erf = [](size_t n, double a, double b, double steepness, const std::vector<double>& expected) {
		const alfi::dist::Erf<double> dist = alfi::dist::Erf<double>(n, a, b, steepness);
		EXPECT_EQ(dist.type(), alfi::dist::Type::ERF);
		EXPECT_EQ(dist.n(), n);
		EXPECT_EQ(dist.a(), a);
		EXPECT_EQ(dist.b(), b);
		EXPECT_EQ(dist.steepness(), steepness);
		EXPECT_EQ(dist.n(), expected.size());
		for (size_t index = 0; index < expected.size(); index++) {
			EXPECT_DOUBLE_EQ(dist[index], expected[index]);
		}
	};

	test_case_erf(0, 0, 1, 1, {});
	test_case_erf(1, 0, 1, 1, {0.5});
}

TEST(DistributionsTest, ErfStretchedDistribution) {

	static const auto test_case_erf_stretched = [](size_t n, double a, double b, double steepness, const std::vector<double>& expected) {
		const alfi::dist::ErfStretched<double> dist = alfi::dist::ErfStretched<double>(n, a, b, steepness);
		EXPECT_EQ(dist.type(), alfi::dist::Type::ERF_STRETCHED);
		EXPECT_EQ(dist.n(), n);
		EXPECT_EQ(dist.a(), a);
		EXPECT_EQ(dist.b(), b);
		EXPECT_EQ(dist.steepness(), steepness);
		EXPECT_LE(dist.stretched_a(), a);
		EXPECT_GE(dist.stretched_b(), b);
		EXPECT_EQ(dist.n(), expected.size());
		for (size_t index = 0; index < expected.size(); index++) {
			EXPECT_DOUBLE_EQ(dist[index], expected[index]);
		}
	};

	test_case_erf_stretched(0, 0, 1, 1, {});
	test_case_erf_stretched(1, 0, 1, 1, {0.5});
	test_case_erf_stretched(2, 0, 1, 1, {0, 1});
	test_case_erf_stretched(3, 0, 1, 1, {0, 0.5, 1});
}