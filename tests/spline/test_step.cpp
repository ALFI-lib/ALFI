#include <ALFI/spline/step.h>

#include "../test_utils.h"

const auto test_data_path = TEST_DATA_DIR "/spline/step.toml";

const auto test_data = toml::parse_file(test_data_path);

void test_step_spline(double epsilon) {
	const auto& test_cases = test_data["test_cases"].ref<toml::array>();

	test_cases.for_each([&](const toml::table& test_case) {
		const auto& type = test_case["type"].ref<std::string>();
		const auto& X = to_vector<double>(test_case["X"].ref<toml::array>());
		const auto& Y = to_vector<double>(test_case["Y"].ref<toml::array>());
		const auto& xx = to_vector<double>(test_case["xx"].ref<toml::array>());
		const auto& yy = to_vector<double>(test_case["yy"].ref<toml::array>());

		const auto t = [&]() -> alfi::spline::StepSpline<>::Type {
			if (type == "left") {
				return alfi::spline::StepSpline<>::Types::Left{};
			} else if (type == "middle") {
				return alfi::spline::StepSpline<>::Types::Middle{};
			} else if (type == "right") {
				return alfi::spline::StepSpline<>::Types::Right{};
			} else {
				throw std::runtime_error{"Unexpected type:" + type};
			}
		}();

		expect_eq(alfi::spline::StepSpline<>(X, Y, t).eval(xx), yy, epsilon);
	});
}

TEST(StepSplineTest, TestData) {
	test_step_spline(1e-15);
}

TEST(StepSplineTest, General) {
	const std::vector<double> X = {0, 1, 2, 3};
	const std::vector<double> Y = {5, 2, 10, 4};
	const auto left_spline = alfi::spline::StepSpline(X, Y, alfi::spline::StepSpline<>::Types::Left{});
	const auto middle_spline = alfi::spline::StepSpline(X, Y, alfi::spline::StepSpline<>::Types::Middle{});
	const auto right_spline = alfi::spline::StepSpline(X, Y, alfi::spline::StepSpline<>::Types::Right{});
	const auto fraction_spline = alfi::spline::StepSpline(X, Y, alfi::spline::StepSpline<>::Types::Fraction{0.25});
	const auto proportion_spline = alfi::spline::StepSpline(X, Y, alfi::spline::StepSpline<>::Types::Proportion{3, 7});
	const auto epsilon = 0.0001;
	for (size_t i = 0; i < X.size() - 1; ++i) {
		// Left exact points
		EXPECT_EQ(left_spline(X[i]), Y[i]);
		EXPECT_EQ(left_spline(X[i+1]), Y[i+1]);
		// Left between points
		EXPECT_EQ(left_spline(X[i] + epsilon), Y[i]);
		EXPECT_EQ(left_spline((X[i]+X[i+1])/2), Y[i]);
		EXPECT_EQ(left_spline(X[i+1] - epsilon), Y[i]);

		// Middle exact points
		EXPECT_EQ(middle_spline(X[i]), Y[i]);
		EXPECT_EQ(middle_spline(X[i+1]), Y[i+1]);
		// Middle between points
		EXPECT_EQ(middle_spline(X[i] + epsilon), (Y[i]+Y[i+1])/2);
		EXPECT_EQ(middle_spline((X[i]+X[i+1])/2), (Y[i]+Y[i+1])/2);
		EXPECT_EQ(middle_spline(X[i+1] - epsilon), (Y[i]+Y[i+1])/2);

		// Right exact points
		EXPECT_EQ(right_spline(X[i]), Y[i]);
		EXPECT_EQ(right_spline(X[i+1]), Y[i+1]);
		// Right between points
		EXPECT_EQ(right_spline(X[i] + epsilon), Y[i+1]);
		EXPECT_EQ(right_spline((X[i]+X[i+1])/2), Y[i+1]);
		EXPECT_EQ(right_spline(X[i+1] - epsilon), Y[i+1]);

		// Fraction exact points
		EXPECT_EQ(fraction_spline(X[i]), Y[i]);
		EXPECT_EQ(fraction_spline(X[i+1]), Y[i+1]);
		// Fraction between points
		EXPECT_EQ(fraction_spline(X[i] + epsilon), Y[i] + 0.25*(Y[i+1] - Y[i]));
		EXPECT_EQ(fraction_spline((X[i]+X[i+1])/2), Y[i] + 0.25*(Y[i+1] - Y[i]));
		EXPECT_EQ(fraction_spline(X[i+1] - epsilon), Y[i] + 0.25*(Y[i+1] - Y[i]));

		// Proportion exact points
		EXPECT_EQ(proportion_spline(X[i]), Y[i]);
		EXPECT_EQ(proportion_spline(X[i+1]), Y[i+1]);
		// Proportion between points
		EXPECT_EQ(proportion_spline(X[i] + epsilon), Y[i] + 0.3*(Y[i+1] - Y[i]));
		EXPECT_EQ(proportion_spline((X[i]+X[i+1])/2), Y[i] + 0.3*(Y[i+1] - Y[i]));
		EXPECT_EQ(proportion_spline(X[i+1] - epsilon), Y[i] + 0.3*(Y[i+1] - Y[i]));
	}
	EXPECT_EQ(left_spline(X.front() - 1), Y.front());
	EXPECT_EQ(left_spline(X.back() + 1), Y.back());
	EXPECT_EQ(middle_spline(X.front() - 1), Y.front());
	EXPECT_EQ(middle_spline(X.back() + 1), Y.back());
	EXPECT_EQ(right_spline(X.front() - 1), Y.front());
	EXPECT_EQ(right_spline(X.back() + 1), Y.back());
	EXPECT_EQ(fraction_spline(X.front() - 1), Y.front());
	EXPECT_EQ(fraction_spline(X.back() + 1), Y.back());
	EXPECT_EQ(proportion_spline(X.front() - 1), Y.front());
	EXPECT_EQ(proportion_spline(X.back() + 1), Y.back());
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