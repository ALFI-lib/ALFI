#pragma once

#include <iostream>

#include "config.h"
#include "dist.h"
#include "util/numeric.h"

namespace alfi::misc {
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> barycentric(
			const Container<Number>& X,
			const Container<Number>& Y,
			const Container<Number>& xx,
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

		// barycentric weights
		Container<Number> W(N);
		if (dist_type == dist::Type::UNIFORM) {
			W[0] = 1;
			for (SizeT j = 1; j <= N / 2; ++j) {
				W[j] = -W[j-1] * (static_cast<Number>(N - j)) / static_cast<Number>(j);
			}
			for (SizeT j = N / 2 + 1; j < N; ++j) {
				if (W[j-1] < 0) {
					W[j] = std::abs(W[N-1-j]);
				} else {
					W[j] = -std::abs(W[N-1-j]);
				}
			}
			// TODO: Investigate if normalization is needed.
			// Normalization is needed because the middle weights grow as O(2^n),
			// while the edges grow as O(1). Dividing by 2^(N/2) balances the exponents,
			// reducing numerical instability.
			// const Number norm_factor = std::pow(static_cast<Number>(2), static_cast<Number>(N / 2));
			// for (SizeT j = 0; j < N; ++j) {
			// 	W[j] /= norm_factor;
			// }
		} else if (dist_type == dist::Type::CHEBYSHEV || dist_type == dist::Type::CHEBYSHEV_STRETCHED) {
			for (SizeT j = 0; j < N; ++j) {
				W[j] = (j % 2 == 0 ? 1 : -1) * std::sin(((2 * static_cast<Number>(j) + 1) * M_PI) / (2 * static_cast<Number>(N)));
			}
		} else if (dist_type == dist::Type::CHEBYSHEV_AUGMENTED) {
			W[0] = 1 / static_cast<Number>(2);
			for (SizeT j = 1; j < N - 1; ++j) {
				W[j] = (j % 2 == 0 ? 1 : -1) / ((N - 2) * std::sin(((static_cast<Number>(2*j - 1) * M_PI) / static_cast<Number>(2*N - 4))));
			}
			W[N-1] = (N % 2 == 0 ? -1 : 1) / static_cast<Number>(2);
		} else if (dist_type == dist::Type::CHEBYSHEV_2) {
			for (SizeT j = 0; j < N; ++j) {
				W[j] = (j % 2 == 0 ? 1 : -1);
				if (j == 0 || j == N - 1) {
					W[j] /= 2;
				}
			}
		} else if (dist_type == dist::Type::CHEBYSHEV_3) {
			for (SizeT j = 0; j < N; ++j) {
				W[j] = (j % 2 == 0 ? 1 : -1) * std::cos((static_cast<Number>(2*j) * M_PI) / static_cast<Number>(2*N - 1) / 2);
				if (j == 0) {
					W[j] /= 2;
				}
			}
		} else if (dist_type == dist::Type::CHEBYSHEV_4) {
			for (SizeT j = 0; j < N; ++j) {
				W[j] = (j % 2 == 0 ? 1 : -1) * std::sin((static_cast<Number>(2*j + 1) * M_PI) / static_cast<Number>(2*N - 1) / 2);
				if (j == N - 1) {
					W[j] /= 2;
				}
			}
		} else {
			// The factor scales the weights to prevent excessive growth or shrinkage, reducing numerical instability.
			// TODO: Investigate if normalization is needed and refine the scaling factor.
			// const Number norm_factor = (X[N-1] - X[0]) / static_cast<Number>(2);
			for (SizeT j = 0; j < N; ++j) {
				W[j] = 1;
				for (SizeT k = 0; k < N; ++k) {
					if (j != k) {
						// c[i] *= norm_factor / (X[j] - X[k]); // with normalization
						W[j] /= (X[j] - X[k]);
					}
				}
			}
		}

		const auto nn = xx.size();

		Container<Number> result(nn);

#if defined(_OPENMP) && !defined(ALFI_DISABLE_OPENMP)
#pragma omp parallel for
#endif
		for (SizeT k = 0; k < nn; ++k) {
			Number numerator = 0, denominator = 0;
			SizeT exact_idx = N;

			for (SizeT i = 0; i < N; ++i) {
				if (util::numeric::are_equal(xx[k], X[i], epsilon)) {
					exact_idx = i;
					break;
				}
				const auto x_diff = xx[k] - X[i];
				const auto temp = W[i] / x_diff;
				numerator += temp * Y[i];
				denominator += temp;
			}

			if (exact_idx != N) {
				result[k] = Y[exact_idx];
			} else {
				result[k] = numerator / denominator;
			}
		}

		return result;
	}
}