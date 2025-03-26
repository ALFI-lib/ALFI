#include <ALFI/misc.h>

#include "../test_utils.h"

void test_barycentric(const toml::parse_result& test_data, double epsilon) {
	const auto& test_cases = test_data["test_cases"].ref<toml::array>();

	test_cases.for_each([&](const toml::table& test_case) {
		const auto& X = to_vector<double>(test_case["X"].ref<toml::array>());
		const auto& Y = to_vector<double>(test_case["Y"].ref<toml::array>());
		const auto& xx = to_vector<double>(test_case["xx"].ref<toml::array>());
		const auto& dist = test_case["dist"].ref<std::string>();
		const auto& yy = to_vector<double>(test_case["yy"].ref<toml::array>());

		const auto& dist_type
			= dist == "uniform" ?
				alfi::dist::Type::UNIFORM
			: dist == "chebyshev" ?
				alfi::dist::Type::CHEBYSHEV
			: dist == "chebyshev_2" ?
				alfi::dist::Type::CHEBYSHEV_2
			: alfi::dist::Type::GENERAL;

		expect_eq(alfi::misc::barycentric(X, Y, xx, dist_type), yy, epsilon);
	});
}

// barycentric formula can also be tested on the polynomial data
TEST(BarycentricTest, PolynomialData) {
	const auto polynomial_test_data_path = TEST_DATA_DIR "/poly/poly.toml";
	const auto polynomial_test_data = toml::parse_file(polynomial_test_data_path);
	test_barycentric(polynomial_test_data, 1e-12);
}

TEST(BarycentricTest, BarycentricData) {
	const auto barycentric_test_data_path = TEST_DATA_DIR "/misc/barycentric.toml";
	const auto barycentric_test_data = toml::parse_file(barycentric_test_data_path);
	test_barycentric(barycentric_test_data, 1e-13);
}