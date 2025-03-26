#include <ALFI/poly.h>

#include "../test_utils.h"

const auto test_data_path = TEST_DATA_DIR "/poly/poly.toml";

const auto test_data = toml::parse_file(test_data_path);

void test_polynomial(const auto& function_coeffs, const auto& function_vals,
					 double epsilon_coeffs, double epsilon_values, double epsilon_vals) {
	const auto& test_cases = test_data["test_cases"].ref<toml::array>();

	test_cases.for_each([&](const toml::table& test_case) {
		const auto& X = to_vector<double>(test_case["X"].ref<toml::array>());
		const auto& Y = to_vector<double>(test_case["Y"].ref<toml::array>());
		const auto& coeffs = to_vector<double>(test_case["coeffs"].ref<toml::array>());
		const auto& xx = to_vector<double>(test_case["xx"].ref<toml::array>());
		const auto& yy = to_vector<double>(test_case["yy"].ref<toml::array>());

		const auto P = function_coeffs(X, Y);
		expect_eq(P, coeffs, epsilon_coeffs);
		const auto values = alfi::poly::val(P, xx);
		expect_eq(values, yy, epsilon_values);
		const auto vals = function_vals(X, Y, xx);
		expect_eq(vals, yy, epsilon_vals);
	});
}

// need this workaround because of the fourth parameter of `alfi::poly::imp_lagrange_vals`
static const auto& imp_lagrange_vals = [](const auto& X, const auto& Y, const auto& xx) {
	return alfi::poly::imp_lagrange_vals<>(X, Y, xx);
};

TEST(PolynomialsTest, Lagrange) {
	test_polynomial(alfi::poly::lagrange<>, alfi::poly::lagrange_vals<>, 1e-7, 1e-1, 1e-11);
}

TEST(PolynomialsTest, ImprovedLagrange) {
	test_polynomial(alfi::poly::imp_lagrange<>, imp_lagrange_vals, 1e-4, 1e-4, 1e-11);
}

TEST(PolynomialsTest, Newton) {
	test_polynomial(alfi::poly::newton<>, alfi::poly::newton_vals<>, 1e-8, 1e-4, 1e-5);
}