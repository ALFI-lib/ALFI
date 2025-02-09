#pragma once

#include <iostream>

#include "../config.h"
#include "../util/linalg.h"

namespace alfi::util::spline {
	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
		Container<Number> simple_spline(const Container<Number>& X, const Container<Number>& Y, SizeT degree) {
		if (X.size() != Y.size()) {
			std::cerr << "Error in function " << __FUNCTION__
					  << ": Vectors X (of size " << X.size()
					  << ") and Y (of size " << Y.size()
					  << ") are not the same size. Returning an empty array..." << std::endl;
			return {};
		}

		const auto n = X.size();
		const auto m = degree + 1;

		if (n > m) {
			std::cerr << "Error in function " << __FUNCTION__
					  << ": Number of points (" << n
					  << ") is bigger than degree + 1 (" << m
					  << "). The spline is not simple. Returning an empty array..." << std::endl;
			return {};
		}

		if (n == 0) {
			return {};
		}

		if (n == 1) {
			Container<Number> coeffs(m, 0);
			coeffs[m-1] = Y[0];
			return coeffs;
		}

		if (n == 2) {
			Container<Number> coeffs(m, 0);
			coeffs[m-2] = (Y[1] - Y[0]) / (X[1] - X[0]);
			coeffs[m-1] = Y[0];
			return coeffs;
		}

		if (n == 3) {
			const auto h1 = X[1] - X[0], h2 = X[2] - X[1];
			const auto d1 = (Y[1] - Y[0]) / h1, d2 = (Y[2] - Y[1]) / h2;

			Container<Number> coeffs(2*m, 0);

			coeffs[m-3] = (d2 - d1) / (h1 + h2);
			coeffs[m-2] = d1 - h1 * coeffs[m-3];
			coeffs[m-1] = Y[0];
			coeffs[2*m-3] = coeffs[m-3];
			coeffs[2*m-2] = d1 + h1 * coeffs[m-3];
			coeffs[2*m-1] = Y[1];

			return coeffs;
		}

		if (n == 4) {
			const auto h1 = X[1]-X[0], h2 = X[2]-X[1];
			const auto h12 = X[2]-X[0], h13 = X[3]-X[0];
			const auto d1 = Y[1]-Y[0], d12 = Y[2]-Y[0], d13 = Y[3]-Y[0];

			const auto abc = util::linalg::lup_solve(
				{{std::pow(h1, 3), std::pow(h1, 2), h1},
					{std::pow(h12, 3), std::pow(h12, 2), h12},
					{std::pow(h13, 3), std::pow(h13, 2), h13}},
				{d1, d12, d13});
			if (abc.empty()) {
				std::cerr << "Error in function " << __FUNCTION__ << ": Could not compute. Returning an empty array..." << std::endl;
				return {};
			}

			const auto a = abc[0], b = abc[1], c = abc[2];

			Container<Number> coeffs(3*m, 0);

			coeffs[m-4] = a;
			coeffs[m-3] = b;
			coeffs[m-2] = c;
			coeffs[m-1] = Y[0];
			coeffs[2*m-4] = a;
			coeffs[2*m-3] = 3 * a * h1 + b;
			coeffs[2*m-2] = 3 * a * std::pow(h1, 2) + 2 * b * h1 + c;
			coeffs[2*m-1] = Y[1];
			coeffs[3*m-4] = a;
			coeffs[3*m-3] = 3 * a * h2 + coeffs[2*m-3];
			coeffs[3*m-2] = 3 * a * std::pow(h2, 2) + 2 * coeffs[2*m-3] * h2 + coeffs[2*m-2];
			coeffs[3*m-1] = Y[2];

			return coeffs;
		}

		std::cerr << "Error in function " << __FUNCTION__
				  << ": Number of points (" << n
				  << ") is too big (bigger than 4). Returning an empty array..." << std::endl;
		return {};
	}
}