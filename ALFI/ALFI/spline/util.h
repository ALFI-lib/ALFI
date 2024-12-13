#pragma once

#include <iostream>

#include "../config.h"

namespace alfi::spline::util {
	template <typename Number = DefaultNumber, typename Container = DefaultContainer<Number>>
		Container simple_spline(const Container& X, const Container& Y, size_t degree) {
		if (X.size() != Y.size()) {
			std::cerr << "Error in function " << __FUNCTION__
					  << ": Vectors X (of size " << X.size()
					  << ") and Y (of size " << Y.size()
					  << ") are not the same size. Returning an empty array..." << std::endl;
			return {};
		}

		const size_t n = X.size();
		const size_t m = degree + 1;

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
			Container coeffs(m, 0);
			coeffs[m-1] = Y[0];
			return coeffs;
		}

		if (n == 2) {
			Container coeffs(m, 0);
			coeffs[m-2] = (Y[1] - Y[0]) / (X[1] - X[0]);
			coeffs[m-1] = Y[0];
			return coeffs;
		}

		if (n == 3) {
			const Number h1 = X[1] - X[0], h2 = X[2] - X[1];
			const Number d1 = (Y[1] - Y[0]) / h1, d2 = (Y[2] - Y[1]) / h2;

			Container coeffs(2*m, 0);

			coeffs[m-3] = (d2 - d1) / (h1 + h2);
			coeffs[m-2] = d1 - h1 * coeffs[m-3];
			coeffs[m-1] = Y[0];
			coeffs[2*m-3] = coeffs[m-3];
			coeffs[2*m-2] = d1 + h1 * coeffs[m-3];
			coeffs[2*m-1] = Y[1];

			return coeffs;
		}

		std::cerr << "Error in function " << __FUNCTION__
				  << ": Number of points (" << n
				  << ") is too big (bigger than 3). Returning an empty array..." << std::endl;
		return {};
	}
}