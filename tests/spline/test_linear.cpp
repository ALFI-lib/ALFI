#include <ALFI/spline/linear.h>

#include "../test_utils.h"

const auto test_data_path = TEST_DATA_DIR "/spline/linear.toml";

const auto test_data = toml::parse_file(test_data_path);

void test_linear_spline(double epsilon_coeffs, double epsilon_values) {
	const auto& test_cases = test_data["test_cases"].ref<toml::array>();

	test_cases.for_each([&](const toml::table& test_case) {
		const auto& X = to_vector<double>(test_case["X"].ref<toml::array>());
		const auto& Y = to_vector<double>(test_case["Y"].ref<toml::array>());
		const auto& coeffs = to_vector<double>(test_case["coeffs"].ref<toml::array>());
		const auto& xx = to_vector<double>(test_case["xx"].ref<toml::array>());
		const auto& yy = to_vector<double>(test_case["yy"].ref<toml::array>());

		const auto spline = alfi::spline::LinearSpline<>(X, Y);
		expect_eq(spline.coeffs(), coeffs, epsilon_coeffs);
		const auto values = spline.eval(xx);
		expect_eq(values, yy, epsilon_values);
	});
}

TEST(LinearSplineTest, TestData) {
	test_linear_spline(1e-14, 1e-14);
}

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