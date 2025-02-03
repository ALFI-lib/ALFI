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

	// WRONG CODE
	// Tried to implement this for "Not-a-knot" spline, because it had the third coefficient in the first equation,
	// but decided that it would be better to derive and equation with only two coefficients.
	// See https://github.com/ALFI-lib/alfi-lib.github.io/pull/8
	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer, typename I1, typename I2, typename I3, typename I4>
	Container<Number> tridiag_solve(
			Number a11, Number a12, Number a13, Number r1,
			Number ann2, Number ann1, Number ann, Number rn,
			I1& lower_begin, const I1& lower_end,
			I2& diag_begin, const I2& diag_end,
			I3& upper_begin, const I3& upper_end,
			I4& right_begin, const I4& right_end
	) {
		const auto n = right_end - right_begin + 2;
		assert(n == (lower_end - lower_begin + 2));
		assert(n == (diag_end - diag_begin + 2));
		assert(n == (upper_end - upper_begin + 2));
		auto lower = lower_begin;
		auto diag = diag_begin;
		auto upper = upper_begin;
		auto right = right_begin;
		const auto m_1 = *lower / a11;
		*diag -= m_1 * a12;
		*upper -= m_1 * a13;
		*right -= m_1 * r1;
		++lower, ++diag, ++upper, ++right;
		for (; lower < lower_end; ++lower, ++diag, ++upper, ++right) {
			const auto m = *lower / *(diag - 1);
			*diag -= m * *(upper - 1);
			*right -= m * *(right - 1);
		}
		if (ann2 != 0) {
			const auto m_n = ann2 / *(lower - 1);
			ann1 -= m_n * *(diag - 1);
			ann -= m_n * *(upper - 1);
			rn -= m_n * *(right - 1);
		}
		const auto m_n = ann1 / *(diag - 1);
		ann -= m_n * *(upper - 1);
		rn -= m_n * *(right - 1);
		Container<Number> X(n);
		X[n-1] = rn / ann;
		for (auto i = n - 2; i > 0; --i) {
			X[i] = (*--right - *--upper * X[i+1]) / *--diag;
		}
		X[0] = (r1 - a12 * X[1] - a13 * X[2]) / a11;
		return X;
	}
}