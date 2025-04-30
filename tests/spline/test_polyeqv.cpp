#include <ALFI/spline/polyeqv.h>

#include "../test_utils.h"

const auto test_data_path = TEST_DATA_DIR "/poly/poly.toml";

const auto test_data = toml::parse_file(test_data_path);

void test_polyeqv_spline(double epsilon_accuracy, double epsilon_speed) {
	const auto& test_cases = test_data["test_cases"].ref<toml::array>();

	test_cases.for_each([&](const toml::table& test_case) {
		const auto& X = to_vector<double>(test_case["X"].ref<toml::array>());
		const auto& Y = to_vector<double>(test_case["Y"].ref<toml::array>());
		const auto& xx = to_vector<double>(test_case["xx"].ref<toml::array>());
		const auto& yy = to_vector<double>(test_case["yy"].ref<toml::array>());

		expect_eq(alfi::spline::PolyEqvSpline<>(X, Y, alfi::spline::PolyEqvSpline<>::OptimizationType::ACCURACY).eval(xx), yy, epsilon_accuracy);
		expect_eq(alfi::spline::PolyEqvSpline<>(X, Y, alfi::spline::PolyEqvSpline<>::OptimizationType::SPEED).eval(xx), yy, epsilon_speed);
	});
}

TEST(PolyEqvSplineTest, PolynomialData) {
	test_polyeqv_spline(1e-10, 1e-2);
}

TEST(PolyEqvSplineTest, General) {
	const std::vector<double> X = {0, 1, 2, 3};
	const std::vector<double> Y = {5, 2, 10, 4};
	const auto spline = alfi::spline::PolyEqvSpline<>(X, Y);
	for (size_t i = 0; i < X.size() - 1; ++i) {
		EXPECT_EQ(spline(X[i]), Y[i]);
	}
	EXPECT_NEAR(spline(X[X.size()-1]), Y[Y.size()-1], 1e-14);
}

TEST(PolyEqvSplineTest, Moving) {
	std::vector<double> X = {0, 1, 2, 3};
	const std::vector<double> Y = {5, 2, 10, 4};
	auto spline = alfi::spline::PolyEqvSpline<>(std::move(X), Y);
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