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

	/**
		@brief Solves a tridiagonal system of linear equations.

		This function solves a system of linear equations \f(AX = B\f), where \f(A\f) is a tridiagonal matrix.
		The input vectors are modified in place.

		Unlike \ref tridiag_solve, this algorithm is generally unstable.\n
		The sufficient conditions for this algorithm to be stable are
		(see this article: https://stsynkov.math.ncsu.edu/book_sample_material/Sections_5.4-5.5.pdf):
			- \f(\abs{diag[0]} \ge \abs{upper[0]}\f),
			- \f(\abs{diag[k]} \ge \abs{lower[k]} + \abs{upper[k]}, k = 1, 2, ..., n - 2\f),
			- \f(\abs{diag[n-1]} \ge \abs{lower[n-1]}\f),
			- one of \f(n\f) inequalities is strict, i.e., \f(>\f) rather than \f(\ge\f).

		The first element of `lower` and the last element of `upper` are ignored.

		@param lower the subdiagonal elements of the tridiagonal matrix (first element is ignored)
		@param diag the diagonal elements of the tridiagonal matrix
		@param upper the superdiagonal elements of the tridiagonal matrix (last element is ignored)
		@param right the right-hand side vector of the system
		@return the solution vector
	 */
	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> tridiag_solve_unstable(const auto& lower, auto&& diag, const auto& upper, auto&& right) {
		const auto n = right.size();
		assert(n == lower.size());
		assert(n == diag.size());
		assert(n == upper.size());
		if (n == 0) {
			return {};
		}
		for (SizeT i = 0; i < n - 1; ++i) {
			const auto m = lower[i+1] / diag[i];
			diag[i+1] -= m * upper[i];
			right[i+1] -= m * right[i];
		}
		Container<Number> X(n);
		X[n-1] = right[n-1] / diag[n-1];
		for (SizeT iter = 2; iter <= n; ++iter) {
			const auto i = n - iter;
			X[i] = (right[i] - upper[i] * X[i+1]) / diag[i];
		}
		return X;
	}

	/**
		@brief Solves a tridiagonal system of linear equations.

		This function solves a system of linear equations \f(AX = B\f), where \f(A\f) is a tridiagonal matrix.
		The input vectors are modified in place.

		Unlike \ref tridiag_solve_unstable, this algorithm is stable.

		The first element of `lower` and the last element of `upper` are ignored.

		@note Based on the https://github.com/snsinfu/cxx-spline/blob/625d4d325cb2/include/spline.hpp#L40-L105 code
			  from https://github.com/snsinfu, which was licensed under the BSL-1.0 (Boost Software License 1.0).

		@param lower the subdiagonal elements of the tridiagonal matrix (first element is ignored)
		@param diag the diagonal elements of the tridiagonal matrix
		@param upper the superdiagonal elements of the tridiagonal matrix (last element is ignored)
		@param right the right-hand side vector of the system
		@return the solution vector
	 */
	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> tridiag_solve(auto&& lower, auto&& diag, auto&& upper, auto&& right) {
		const auto n = right.size();
		assert(n == lower.size());
		assert(n == diag.size());
		assert(n == upper.size());
		if (n == 0) {
			return {};
		}
		if (n == 1) {
			return {right[0]/diag[0]};
		}
		for (SizeT i = 0; i < n - 1; ++i) {
			if (std::abs(diag[i]) >= std::abs(lower[i+1])) {
				const auto m = lower[i+1] / diag[i];
				diag[i+1] -= m * upper[i];
				right[i+1] -= m * right[i];
				lower[i+1] = 0;
			} else {
				// Swap rows i and (i+1).
				// Eliminate the lower[i+1] element by reducing row (i+1) by row i.
				// Use lower[i+1] as a buffer for the non-tridiagonal element (i,i+2) (hack, used below).
				const auto m = diag[i] / lower[i+1];
				diag[i] = lower[i+1];
				lower[i+1] = upper[i+1];
				upper[i+1] *= -m;
				std::swap(upper[i], diag[i+1]);
				diag[i+1] -= m * upper[i];
				std::swap(right[i], right[i+1]);
				right[i+1] -= m * right[i];
			}
		}
		Container<Number> X(n);
		X[n-1] = right[n-1] / diag[n-1];
		X[n-2] = (right[n-2] - upper[n-2] * X[n-1]) / diag[n-2];
		for (SizeT iter = 3; iter <= n; ++iter) {
			const auto i = n - iter;
			X[i] = (right[i] - upper[i] * X[i+1] - lower[i+1] * X[i+2]) / diag[i];
		}
		return X;
	}
}