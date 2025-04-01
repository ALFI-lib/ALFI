#pragma once

#include <cassert>
#include <cmath>
#include <iostream>

#include "../config.h"

namespace alfi::util::linalg {
	/**
		@brief Solves a system of linear equations using LUP decomposition.

		This function solves the system of linear equations \f(AX = B\f), where \f(A\f) is a square matrix,
		by performing LUP decomposition and forward-backward substitution.\n
		The input matrix \f(A\f) is decomposed into lower \f(L\f) and upper \f(U\f) triangular matrices,
		such that \f(PA = LU\f), transforming the equation into \f(PLUX = B\f) or \f(LUX = PB\f).

		The solution is found in two stages:
		1. The first stage solves the system \f(LY = PB\f) for the intermediate vector \f(Y\f) using forward substitution.
		2. The second stage solves the system \f(UX = Y\f) for the final solution \f(X\f) using backward substitution.

		If the pivot element in a given column is too small (less than the specified \f(epsilon\f)),
		the function returns an empty container, indicating that the matrix is degenerate,
		and prints the corresponding message to `std::cerr`.

		The input matrix \f(A\f) and vector \f(B\f) are modified in place.

		@param A the matrix (2D container of size \f(n \times n\f))
		@param B the right-hand side vector (1D container of size \f(n\f))
		@param epsilon a small threshold value to detect degeneracy (default is machine epsilon)
		@return the solution vector or an empty container if the matrix is degenerate
	 */
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> lup_solve(Container<Container<Number>>&& A, Container<Number>&& B, Number epsilon = std::numeric_limits<Number>::epsilon()) {
		const auto n = B.size();
		assert(n == A.size());
		for (SizeT i = 0; i < n; ++i) {
			assert(n == A[i].size());
			Number max_a = std::abs(A[i][i]);
			SizeT i_max = i;
			for (SizeT k = i + 1; k < n; ++k) {
				const auto cur = std::abs(A[k][i]);
				if (cur > max_a) {
					max_a = cur;
					i_max = k;
				}
			}
			if (max_a <= epsilon) {
				std::cerr << "Error in function " << __FUNCTION__ << ": Matrix A is degenerate. Returning an empty array..." << std::endl;
				return {};
			}
			if (i_max != i) {
				std::swap(A[i], A[i_max]);
				std::swap(B[i], B[i_max]);
			}
			for (SizeT j = i + 1; j < n; ++j) {
				A[j][i] /= A[i][i];
				for (SizeT k = i + 1; k < n; ++k) {
					A[j][k] -= A[j][i] * A[i][k];
				}
			}
		}
		Container<Number> X(n);
		for (SizeT i = 0; i < n; ++i) {
			X[i] = B[i];
			for (SizeT k = 0; k < i; ++k) {
				X[i] -= A[i][k] * X[k];
			}
		}
		for (SizeT iter = 0; iter < n; ++iter) {
			const auto i = n - 1 - iter;
			for (SizeT k = i + 1; k < n; ++k) {
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
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> tridiag_solve_unstable(
			const Container<Number>& lower,
			Container<Number>&& diag,
			const Container<Number>& upper,
			Container<Number>&& right
	) {
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
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> tridiag_solve(
			Container<Number>&& lower,
			Container<Number>&& diag,
			Container<Number>&& upper,
			Container<Number>&& right
	) {
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