#pragma once

#include <iostream>

#include "config.h"
#include "dist.h"
#include "util/numeric.h"

namespace alfi::misc {

	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> barycentric(
			const auto& X,
			const auto& Y,
			const auto& xx,
			dist::Type dist_type = dist::Type::GENERAL,
			Number epsilon = std::numeric_limits<Number>::epsilon()
	) {
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

		const auto N = X.size();

		Container<Number> c(N);
		if (dist_type == dist::Type::UNIFORM) {
			c[0] = 1;
			for (SizeT j = 1; j <= N / 2; ++j) {
				c[j] = -c[j-1] * (static_cast<Number>(N - j)) / static_cast<Number>(j);
			}
			for (SizeT j = N / 2 + 1; j < N; ++j) {
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
			// for (SizeT j = 0; j < N; ++j) {
			// 	c[j] /= norm_factor;
			// }
		} else if (dist_type == dist::Type::CHEBYSHEV) {
			for (SizeT j = 0; j < N; ++j) {
				c[j] = (j % 2 == 0 ? 1 : -1) * sin(((2 * static_cast<Number>(j) + 1) * M_PI) / (2 * static_cast<Number>(N)));
			}
		} else if (dist_type == dist::Type::CIRCLE_PROJECTION) {
			for (SizeT j = 0; j < N; ++j) {
				c[j] = (j % 2 == 0 ? 1 : -1) * (j == 0 || j == N - 1 ? static_cast<Number>(1)/static_cast<Number>(2) : static_cast<Number>(1));
			}
		} else {
			// The factor scales the coefficients to prevent excessive growth or shrinkage, reducing numerical instability.
			// TODO: Investigate if normalization is needed and refine the scaling factor.
			// const Number norm_factor = (X[N-1] - X[0]) / static_cast<Number>(2);
			for (SizeT j = 0; j < N; ++j) {
				c[j] = 1;
				for (SizeT k = 0; k < N; ++k) {
					if (j != k) {
						// c[i] *= norm_factor / (X[j] - X[k]); // with normalization
						c[j] /= (X[j] - X[k]);
					}
				}
			}
		}

		const auto nn = xx.size();

		Container<Number> result(nn);

		const auto NOT_EXACT = N;

		for (SizeT k = 0; k < nn; ++k) {
			Number numerator = 0, denominator = 0;
			SizeT exact = NOT_EXACT;

			for (SizeT i = 0; i < N; ++i) {
				if (util::numeric::are_equal(xx[k], X[i], epsilon)) {
					exact = i;
					break;
				}
				const auto x_diff = xx[k] - X[i];
				const auto temp = c[i] / x_diff;
				numerator += temp * Y[i];
				denominator += temp;
			}

			if (exact != NOT_EXACT) {
				result[k] = Y[exact];
			} else {
				result[k] = numerator / denominator;
			}
		}

		return result;
	}
}