#pragma once

#include "../config.h"

namespace alfi::poly {
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	std::enable_if_t<!traits::has_size<Number>::value, Number>
	val(const Container<Number>& coeffs, const Number& x) {
		Number result = 0;
		for (const Number& c : coeffs) {
			result = result * x + c;
		}
		return result;
	}

	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	std::enable_if_t<traits::has_size<Container<Number>>::value, Container<Number>>
	val(const Container<Number>& coeffs, const Container<Number>& xx) {
		Container<Number> result(xx.size(), 0);
#if defined(_OPENMP) && !defined(ALFI_DISABLE_OPENMP)
#pragma omp parallel for
#endif
		for (SizeT i = 0; i < xx.size(); ++i) {
			for (const Number& c : coeffs) {
				result[i] = result[i] * xx[i] + c;
			}
		}
		return result;
	}
}