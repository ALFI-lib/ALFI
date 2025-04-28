#pragma once

#include <algorithm>
#include <limits>

#include "../config.h"

/**
	@namespace alfi::util::poly
	@brief Namespace with utility functions for polynomial operations.

	This namespace provides a set of functions for manipulating polynomials, such as normalization, multiplication, and division.
	All polynomials are represented as containers of coefficients in descending degree order.
*/
namespace alfi::util::poly {
	/**
		@brief Normalizes a polynomial by removing leading zero coefficients.

		Removes leading coefficients that are considered zero within a given tolerance @p epsilon.
		If all coefficients are close to zero, the last one is preserved.
		If the container is empty, a single zero coefficient is inserted.

		@param p the polynomial to normalize (a container of coefficients in descending degree order)
		@param epsilon the tolerance used to determine whether a coefficient is considered zero (default is machine epsilon)
	 */
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	void normalize(Container<Number>& p, const Number& epsilon = std::numeric_limits<Number>::epsilon()) {
		if (p.empty()) {
			return p.push_back(0);
		}
		auto p_start = std::find_if(p.begin(), p.end(), [&epsilon](Number v) { return std::abs(v) > epsilon; });
		if (p_start == p.end()) {
			--p_start;
		}
		if (p_start > p.begin()) {
			p.erase(p.begin(), p_start);
		}
	}

	/**
		@brief Multiplies two polynomials.

		Given two polynomials \f(p_1\f) and \f(p_2\f) represented as containers of coefficients,
		this function computes their product using the convolution formula:
		\f[
			(p_1 \cdot p_2)[k] = \sum_{i+j=k}{p_1[i] \cdot p_2[j]}
		\f]

		If either polynomial is empty, the function returns an empty container.

		@param p1 the first polynomial
		@param p2 the second polynomial
		@return the product polynomial (either empty or of size `p1.size() + p2.size() - 1`)
	*/
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> mul(const Container<Number>& p1, const Container<Number>& p2) {
		if (p1.empty() || p2.empty()) {
			return {};
		}
		Container<Number> result(p1.size() + p2.size() - 1);
		for (SizeT i = 0; i < p1.size(); ++i) {
			for (SizeT j = 0; j < p2.size(); ++j) {
				result[i+j] += p1[i] * p2[j];
			}
		}
		return result;
	}

	/**
		@brief Divides one polynomial by another.

		This function divides the @p dividend polynomial \f(A\f) by the @p divisor polynomial \f(B\f),
		and returns the @p quotient \f(Q\f) and the @p remainder \f(R\f) such that:
		\f[
			A(x) = B(X) \cdot Q(x) + R(x)
		\f]
		with either \f(R\f) being effectively zero or the degree of \f(R\f) being lower than the degree of \f(B\f).

		The division is performed using a tolerance @p epsilon to determine when a coefficient is considered zero.

		If the divisor is effectively zero or if the dividend has a lower degree than the divisor,
		the function returns an empty quotient and the dividend as the remainder.

		@param dividend the polynomial to be divided
		@param divisor the polynomial to divide by
		@param epsilon the tolerance used to determine whether a coefficient is considered zero (default is machine epsilon)
		@return a pair `{quotient, remainder}`
	 */
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	std::pair<Container<Number>,Container<Number>> div(const Container<Number>& dividend, const Container<Number>& divisor, const Number& epsilon = std::numeric_limits<Number>::epsilon()) {
		const auto divisor_start = std::find_if(divisor.begin(), divisor.end(), [&epsilon](Number v) { return std::abs(v) > epsilon; });

		if (divisor_start == divisor.end() || dividend.size() < divisor.size()) {
			return {{}, dividend};
		}

		const auto divisor_start_idx = divisor_start - divisor.begin();

		const auto n = dividend.size();
		const auto m = divisor.size() - divisor_start_idx;

		Container<Number> quotient(n - m + 1, 0);
		Container<Number> remainder = dividend;

		for (SizeT i = 0; i <= n - m; ++i) {
			const Number factor = remainder[i] / divisor[divisor_start_idx];
			quotient[i] = factor;
			for (SizeT j = 0; j < m; ++j) {
				remainder[i+j] -= factor * divisor[divisor_start_idx+j];
			}
		}

		remainder.erase(remainder.begin(), remainder.end() - m + 1);
		return {quotient, remainder};
	}
}