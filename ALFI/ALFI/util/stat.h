#pragma once

#include "../config.h"

#include <cmath>

namespace alfi::util::stat {
	/**
		@brief Calculates the arithmetic mean of the elements in a container.

		The arithmetic mean is calculated as:
		\f[
			mean = \frac{1}{n} \sum_{i=1}^{n} x_i
		\f]
		where \f(n\f) is the number of elements in the container, and \f(x_i\f) are the elements.

		If the container is empty, the function returns 0.

		@param container the container holding the elements
		@return arithmetic mean of the elements
	*/
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Number mean(const Container<Number>& container) {
		if (container.empty()) {
			return 0;
		}
		Number sum = 0;
		for (const auto& value : container) {
			sum += value;
		}
		return sum / container.size();
	}

	/**
		@brief Calculates the arithmetic mean of the absolute values of the elements in a container.

		The mean absolute value is calculated as:
		\f[
			mean\_abs = \frac{1}{n} \sum_{i=1}^{n} |x_i|
		\f]
		where \f(n\f) is the number of elements in the container, and \f(x_i\f) are the elements.

		If the container is empty, the function returns 0.

		@param container the container holding the elements
		@return arithmetic mean of the absolute values of the elements
	*/
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Number mean_abs(const Container<Number>& container) {
		if (container.empty()) {
			return 0;
		}
		Number sum = 0;
		for (const auto& value : container) {
			sum += std::abs(value);
		}
		return sum / container.size();
	}

	/**
		@brief Calculates the arithmetic mean the squares of the elements in a container.

		The mean square is calculated as:
		\f[
			mean\_square = \frac{1}{n} \sum_{i=1}^{n} x_i^2
		\f]
		where \f(n\f) is the number of elements in the container, and \f(x_i\f) are the elements.

		If the container is empty, the function returns 0.

		@param container the container holding the elements
		@return arithmetic mean the squares of the elements
	*/
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Number mean_square(const Container<Number>& container) {
		if (container.empty()) {
			return 0;
		}
		Number sum = 0;
		for (const auto& value : container) {
			sum += value * value;
		}
		return sum / container.size();
	}

	/**
		@brief Calculates the root mean square (RMS) of the elements in a container.

		The root mean square is calculated as:
		\f[
			rms = \sqrt{\frac{1}{n} \sum_{i=1}^{n} x_i^2}
		\f]
		where \f(n\f) is the number of elements in the container, and \f(x_i\f) are the elements.

		If the container is empty, the function returns 0.

		@param container the container holding the elements
		@return root mean square of the elements
	*/
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Number rms(const Container<Number>& container) {
		return std::sqrt(mean_square(container));
	}

	/**
		@brief Calculates the variance of the elements in a container.

		The variance is calculated as:
		\f[
			var = \frac{1}{n-1} \sum_{i=1}^{n} (x_i - \mu)^2
		\f]
		for an unbiased estimate, where \f(n\f) is the number of elements, \f(x_i\f) are the elements, and \f(\mu\f) is the mean of the elements.

		For a biased estimate, the formula is:
		\f[
			var = \frac{1}{n} \sum_{i=1}^{n} (x_i - \mu)^2
		\f]

		If the container has fewer than two elements, the function returns 0.

		@param container the container holding the elements
		@param biased If false (default), calculates the unbiased variance (divides by \f(n-1\f)); if true, calculates the biased variance (divides by \f(n\f)).
		@return variance of the elements
	*/
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Number var(const Container<Number>& container, bool biased = false) {
		if (container.size() <= 1) {
			return 0;
		}
		const auto m = mean(container);
		Number sum = 0;
		for (const auto& value : container) {
			sum += (value - m) * (value - m);
		}
		return sum / (biased ? container.size() : container.size() - 1);
	}
}