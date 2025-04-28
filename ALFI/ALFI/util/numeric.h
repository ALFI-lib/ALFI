#pragma once

#include <cmath>

#include "../config.h"

namespace alfi::util::numeric {
	template <typename Number = DefaultNumber>
	bool are_equal(const Number& a, const Number& b, const Number& epsilon = std::numeric_limits<Number>::epsilon()) {
		return std::abs(a - b) <= epsilon || std::abs(a - b) <= std::max(std::abs(a), std::abs(b)) * epsilon;
	}

	/**
		@brief Computes the power of a number using binary exponentiation.

		Calculates \f(x^n\f) in \f(O(\log{n})\f) operations using the binary (exponentiation by squaring) method.

		It supports both signed and unsigned exponent types (@p ExponentType).@n
		If the exponent is negative, the function computes the reciprocal of the positive exponentiation.

		@param x the base
		@param n the exponent
		@return \f(x^n\f)
	 */
	template <typename Number, typename ExponentType>
	Number binpow(Number x, ExponentType n) {
		if constexpr (std::is_signed_v<ExponentType>) {
			if (n < 0) {
				return 1 / binpow(x, -n);
			}
		}
		Number result = 1;
		while (n > 0) {
			if ((n & 1) == 1) {
				result *= x;
			}
			x *= x;
			n >>= 1;
		}
		return result;
	}
}