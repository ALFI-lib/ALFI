#pragma once

#include <iostream>

#include "../config.h"

namespace alfi::util::linalg {
	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> linsolve(Container<Container<Number>>&& A, const Container<Number>& B, Number epsilon = std::numeric_limits<Number>::epsilon()) {
		const auto N = B.size();

		Container<SizeT> P(N);
		for (SizeT i = 0; i < N; ++i) {
			P[i] = i;
		}

		for (SizeT i = 0; i < N; ++i) {
			Number max_a = std::abs(A[i][i]);
			SizeT i_max = i;

			for (SizeT k = i + 1; k < N; ++k) {
				const auto cur = std::abs(A[k][i]);
				if (cur > max_a) {
					max_a = cur;
					i_max = k;
				}
			}

			if (max_a < epsilon) {
				std::cerr << "Error in function " << __FUNCTION__ << ": Matrix A is degenerate. Returning an empty array..." << std::endl;
				return {};
			}

			if (i_max != i) {
				std::swap(P[i], P[i_max]);
				std::swap(A[i], A[i_max]);
			}

			for (SizeT j = i + 1; j < N; ++j) {
				A[j][i] /= A[i][i];
				for (SizeT k = i + 1; k < N; ++k) {
					A[j][k] -= A[j][i] * A[i][k];
				}
			}
		}

		Container<Number> X(N);

		for (SizeT i = 0; i < N; ++i) {
			X[i] = B[P[i]];
			for (SizeT k = 0; k < i; ++k) {
				X[i] -= A[i][k] * X[k];
			}
		}

		for (SizeT iter = 0; iter < N; ++iter) {
			const auto i = N - 1 - iter;
			for (SizeT k = i + 1; k < N; ++k) {
				X[i] -= A[i][k] * X[k];
			}
			X[i] /= A[i][i];
		}

		return X;
	}
}