#include <ALFI/spline/quadratic.h>
#include <ALFI/dist.h>

#include "../test_utils.h"

const auto test_data_path = TEST_DATA_DIR "/spline/quadratic.toml";

const auto test_data = toml::parse_file(test_data_path);

void test_quadratic_spline(double epsilon_coeffs, double epsilon_values) {
	const auto& test_cases = test_data["test_cases"].ref<toml::array>();

	test_cases.for_each([&](const toml::table& test_case) {
		const auto& type = test_case["type"].ref<std::string>();
		const auto& X = to_vector<double>(test_case["X"].ref<toml::array>());
		const auto& Y = to_vector<double>(test_case["Y"].ref<toml::array>());
		const auto& coeffs = to_vector<double>(test_case["coeffs"].ref<toml::array>());
		const auto& xx = to_vector<double>(test_case["xx"].ref<toml::array>());
		const auto& yy = to_vector<double>(test_case["yy"].ref<toml::array>());

		const auto t = [&]() -> alfi::spline::QuadraticSpline<>::Type {
			if (type == "semi-not-a-knot") {
				return alfi::spline::QuadraticSpline<>::Types::SemiNotAKnot{};
			} else if (type == "semi-natural") {
				return alfi::spline::QuadraticSpline<>::Types::SemiNatural{};
			} else {
				throw std::runtime_error{"Unexpected type:" + type};
			}
		}();

		const auto spline = alfi::spline::QuadraticSpline<>(X, Y, t);
		expect_eq(spline.coeffs(), coeffs, epsilon_coeffs);
		const auto values = spline.eval(xx);
		expect_eq(values, yy, epsilon_values);
	});
}

TEST(QuadraticSplineTest, TestData) {
	test_quadratic_spline(1e-14, 1e-14);
}

TEST(QuadraticSplineTest, General) {
	const std::vector<double> X = {0, 1, 2, 3};
	const std::vector<double> Y = {5, 2, 10, 4};
	const auto spline = alfi::spline::QuadraticSpline(X, Y);
	for (size_t i = 0; i < X.size() - 1; ++i) {
		EXPECT_EQ(spline(X[i]), Y[i]);
		EXPECT_EQ(spline(X[i+1]), Y[i+1]);
	}
}

TEST(QuadraticSplineTest, Conditions) {
	const size_t n = 10;
	const double a = -1, b = 1;
	const auto X = alfi::dist::uniform(n, a, b);
	const auto dX = alfi::util::arrays::diff(X);
	const auto Y = alfi::dist::chebyshev_2(n, a, b);
	// natural-start, natural-end, semi-natural
	auto spline1 = alfi::spline::QuadraticSpline<>(X, Y, alfi::spline::QuadraticSpline<>::Types::NaturalStart{});
	EXPECT_NEAR(0, 2*spline1.coeffs()[0], 1e-17);
	auto spline2 = alfi::spline::QuadraticSpline<>(X, Y, alfi::spline::QuadraticSpline<>::Types::NaturalEnd{});
	EXPECT_NEAR(0, 2*spline2.coeffs()[3*8], 1e-17);
	const auto spline3 = alfi::spline::QuadraticSpline<>(X, Y, alfi::spline::QuadraticSpline<>::Types::SemiNatural{});
	ASSERT_EQ(spline1.coeffs().size(), spline3.coeffs().size());
	ASSERT_EQ(spline2.coeffs().size(), spline3.coeffs().size());
	for (size_t i = 0; i < spline3.coeffs().size(); ++i) {
		EXPECT_EQ(spline3.coeffs()[i], (spline1.coeffs()[i] + spline2.coeffs()[i]) / 2);
	}
	// not-a-knot-start, not-a-knot-end, semi-not-a-knot
	spline1.construct(X, Y, alfi::spline::QuadraticSpline<>::Types::NotAKnotStart{});
	EXPECT_EQ(2*spline1.coeffs()[0], 2*spline1.coeffs()[3*1]);
	spline2.construct(X, Y, alfi::spline::QuadraticSpline<>::Types::NotAKnotEnd{});
	EXPECT_EQ(2*spline2.coeffs()[3*7], 2*spline2.coeffs()[3*8]);
	const auto spline4 = alfi::spline::QuadraticSpline<>(X, Y, alfi::spline::QuadraticSpline<>::Types::SemiNotAKnot{});
	ASSERT_EQ(spline1.coeffs().size(), spline4.coeffs().size());
	ASSERT_EQ(spline2.coeffs().size(), spline4.coeffs().size());
	for (size_t i = 0; i < spline4.coeffs().size(); ++i) {
		EXPECT_EQ(spline4.coeffs()[i], (spline1.coeffs()[i] + spline2.coeffs()[i]) / 2);
	}
	// semi-semi
	const auto spline5 = alfi::spline::QuadraticSpline<>(X, Y, alfi::spline::QuadraticSpline<>::Types::SemiSemi{});
	ASSERT_EQ(spline3.coeffs().size(), spline5.coeffs().size());
	ASSERT_EQ(spline4.coeffs().size(), spline5.coeffs().size());
	for (size_t i = 0; i < spline5.coeffs().size(); ++i) {
		EXPECT_NEAR(spline5.coeffs()[i], (spline3.coeffs()[i] + spline3.coeffs()[i]) / 2, 1e-15);
	}
	// clamped
	spline1.construct(X, Y, alfi::spline::QuadraticSpline<>::Types::ClampedStart{10});
	EXPECT_EQ(10, spline1.coeffs()[1]);
	spline2.construct(X, Y, alfi::spline::QuadraticSpline<>::Types::ClampedEnd{10});
	EXPECT_NEAR(10, 2*spline2.coeffs()[3*8]*dX[8] + spline2.coeffs()[3*8+1], 1e-14);
	const auto spline6 = alfi::spline::QuadraticSpline<>(X, Y, alfi::spline::QuadraticSpline<>::Types::SemiClamped{10, 10});
	ASSERT_EQ(spline1.coeffs().size(), spline6.coeffs().size());
	ASSERT_EQ(spline2.coeffs().size(), spline6.coeffs().size());
	for (size_t i = 0; i < spline6.coeffs().size(); ++i) {
		EXPECT_EQ(spline6.coeffs()[i], (spline1.coeffs()[i] + spline2.coeffs()[i]) / 2);
	}
	spline1.construct(X, Y, alfi::spline::QuadraticSpline<>::Types::Clamped{5, 10});
	EXPECT_EQ(10, spline1.coeffs()[3*5+1]);
	EXPECT_NEAR(10, 2*spline1.coeffs()[3*4]*dX[8] + spline1.coeffs()[3*4+1], 1e-15);
	// fixed-second
	spline1.construct(X, Y, alfi::spline::QuadraticSpline<>::Types::FixedSecondStart{10});
	EXPECT_EQ(10, 2*spline1.coeffs()[0]);
	spline2.construct(X, Y, alfi::spline::QuadraticSpline<>::Types::FixedSecondEnd{10});
	EXPECT_EQ(10, 2*spline2.coeffs()[3*8]);
	const auto spline7 = alfi::spline::QuadraticSpline<>(X, Y, alfi::spline::QuadraticSpline<>::Types::SemiFixedSecond{10, 10});
	ASSERT_EQ(spline1.coeffs().size(), spline7.coeffs().size());
	ASSERT_EQ(spline2.coeffs().size(), spline7.coeffs().size());
	for (size_t i = 0; i < spline7.coeffs().size(); ++i) {
		EXPECT_EQ(spline7.coeffs()[i], (spline1.coeffs()[i] + spline2.coeffs()[i]) / 2);
	}
	spline1.construct(X, Y, alfi::spline::QuadraticSpline<>::Types::FixedSecond{5, 10});
	EXPECT_EQ(10, 2*spline1.coeffs()[3*5]);
	// not-a-knot
	spline1.construct(X, Y, alfi::spline::QuadraticSpline<>::Types::NotAKnot{2});
	EXPECT_EQ(2*spline1.coeffs()[3*1], 2*spline1.coeffs()[3*2]);
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