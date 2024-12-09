#pragma once

#include <iostream>

#include "config.h"
#include "dist.h"

namespace alfi::misc {

	template <typename Number = DefaultNumber,typename Container = DefaultContainer<Number>, Number epsilon = std::numeric_limits<Number>::epsilon()>
	Container barycentric(const Container& X, const Container& Y, const Container& xx, dist::Type dist_type = dist::Type::GENERAL) {
		if (X.size() != Y.size()) {
			std::cerr << "Error in function " << __FUNCTION__
					  << ": Vectors X (of size " << X.size()
					  << ") and Y (of size " << Y.size()
					  << ") are not the same size. Returning an empty array..." << std::endl;
			return {};
		}

		if (X.empty()) {
			std::cerr << "Error in function " << __FUNCTION__
					  << ": Vectors X and Y are empty. Cannot interpolate. Returning an empty array..." << std::endl;
			return {};	
		}

		const size_t N = X.size();

		Container c(N);
		if (dist_type == dist::Type::UNIFORM) {
			c[0] = 1;
			for (size_t j = 1; j <= N / 2; ++j) {
				c[j] = -c[j-1] * (static_cast<Number>(N - j)) / static_cast<Number>(j);
			}
			for (size_t j = N / 2 + 1; j < N; ++j) {
				if (c[j-1] < 0) {
					c[j] = std::abs(c[N-1-j]);
				} else {
					c[j] = -std::abs(c[N-1-j]);
				}
			}
			// TODO: Investigate if normalization is needed.
			// Normalization is needed because the middle coefficients grow as O(2^n), 
			// while the edges grow as O(1). Dividing by 2^(N/2) balances the exponents,
			// reducing numerical instability.
			// const Number norm_factor = std::pow(static_cast<Number>(2), static_cast<Number>(N / 2));
			// for (size_t j = 0; j < N; ++j) {
			// 	c[j] /= norm_factor;
			// }
		} else if (dist_type == dist::Type::CHEBYSHEV) {
			for (size_t j = 0; j < N; ++j) {
				c[j] = (j % 2 == 0 ? 1 : -1) * sin(((2 * static_cast<Number>(j) + 1) * M_PI) / (2 * static_cast<Number>(N)));
			}
		} else if (dist_type == dist::Type::CIRCLE_PROJECTION) {
			for (size_t j = 0; j < N; ++j) {
				c[j] = (j % 2 == 0 ? 1 : -1) * (j == 0 || j == N - 1 ? static_cast<Number>(1)/static_cast<Number>(2) : static_cast<Number>(1));
			}
		} else {
			// The factor scales the coefficients to prevent excessive growth or shrinkage, reducing numerical instability.
			// TODO: Investigate if normalization is needed and refine the scaling factor.
			// const Number norm_factor = (X[N-1] - X[0]) / static_cast<Number>(2);
			for (size_t j = 0; j < N; ++j) {
				c[j] = 1;
				for (size_t k = 0; k < N; ++k) {
					if (j != k) {
						// c[i] *= norm_factor / (X[j] - X[k]); // with normalization
						c[j] /= (X[j] - X[k]);
					}
				}
			}
		}

		const size_t nn = xx.size();

		Container numerators(nn, static_cast<Number>(0));
		Container denominators(nn, static_cast<Number>(0));

		const size_t NOT_EXACT = N;
		DefaultContainer<size_t> exact(nn, NOT_EXACT);

		for (size_t k = 0; k < N; ++k) {
			for (size_t i = 0; i < nn; ++i) {
				const Number x_diff = xx[i] - X[k];
				if (std::abs(x_diff) < epsilon) {
					exact[i] = static_cast<Number>(k);
					continue;
				}
				const Number temp = c[k] / x_diff;
				numerators[i] += temp * Y[k];
				denominators[i] += temp;
			}
		}

		Container result(nn, static_cast<Number>(0));
		for (size_t i = 0; i < nn; ++i) {
			if (exact[i] != NOT_EXACT) {
				result[i] = Y[exact[i]];
			} else {
				result[i] = numerators[i] / denominators[i];
			}
		}

		return result;
	}
}