#pragma once

#include "config.h"

#include <cmath>

namespace alfi::dist {
	enum class Type {
		GENERAL,
		UNIFORM,
		CHEBYSHEV,
		CHEBYSHEV_STRETCHED,
		CHEBYSHEV_ELLIPSE,
		CHEBYSHEV_ELLIPSE_STRETCHED,
		CIRCLE_PROJECTION,
		ELLIPSE_PROJECTION,
		SIGMOID,
		SIGMOID_STRETCHED,
		ERF,
		ERF_STRETCHED,
	};

	template <typename Number = DefaultNumber>
	void stretch(auto& points, Number a, Number b) {
		if (points.empty()) {
			return;
		}

		if (points.size() == 1 || points.front() == points.back()) {
			std::fill(points.begin(), points.end(), (a + b) / 2);
			return;
		}

		const Number& min = points.front();
		const Number& max = points.back();
		const Number mid = (a + b) / 2;
		const Number scale = (b - a) / (max - min);

		for (SizeT i = 1; i < points.size() - 1; ++i) {
			points[i] = mid + scale * (points[i] - mid);
		}

		points.front() = a;
		points.back() = b;
	}

	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> uniform(SizeT n, Number a, Number b) {
		if (n == 1)
			return {(a+b)/2};
		Container<Number> points(n);
		for (SizeT i = 0; i < n; ++i) {
			points[i] = a + (b - a) * static_cast<Number>(i) / static_cast<Number>(n - 1);
		}
		return points;
	}

	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> chebyshev(SizeT n, Number a, Number b) {
		if (n == 1)
			return {(a+b)/2};
		Container<Number> points(n);
		for (SizeT i = 0; i < n; ++i) {
			const Number x = std::cos(M_PI * (2 * (static_cast<Number>(n) - static_cast<Number>(i)) - 1) / (2 * static_cast<Number>(n)));
			points[i] = a + (x + 1) * (b - a) / 2;
		}
		return points;
	}

	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> chebyshev_stretched(SizeT n, Number a, Number b) {
		Container<Number> points = chebyshev(n, a, b);
		stretch(points, a, b);
		return points;
	}

	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> chebyshev_ellipse(SizeT n, Number a, Number b, Number ratio) {
		Container<Number> points(n);
		for (SizeT i = 0; i < n / 2; ++i) {
			const Number theta = M_PI * (2 * static_cast<Number>(i) + 1) / (2 * static_cast<Number>(n));
			const Number x = 1 / std::sqrt(1 + std::pow(std::tan(theta) / ratio, 2));
			points[i] = a + (1 - x) * (b - a) / 2;
			points[n-1-i] = a + (1 + x) * (b - a) / 2;
		}
		if (n % 2 == 1)
			points[n/2] = (a+b)/2;
		return points;
	}

	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> chebyshev_ellipse_stretched(SizeT n, Number a, Number b, Number ratio) {
		Container<Number> points = chebyshev_ellipse(n, a, b, ratio);
		stretch(points, a, b);
		return points;
	}

	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> circle_proj(SizeT n, Number a, Number b) {
		if (n == 1)
			return {(a+b)/2};
		Container<Number> points(n);
		for (SizeT i = 0; i < n; ++i) {
			const Number x = 1 - std::cos(M_PI * static_cast<Number>(i) / (static_cast<Number>(n) - 1));
			points[i] = a + (b - a) * x / 2;
		}
		return points;
	}

	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> ellipse_proj(SizeT n, Number a, Number b, Number ratio) {
		Container<Number> points(n);
		for (SizeT i = 0; i < n / 2; ++i) {
			const Number theta = M_PI * static_cast<Number>(i) / (static_cast<Number>(n) - 1);
			const Number x = 1 / std::sqrt(1 + std::pow(std::tan(theta) / ratio, 2));
			points[i] = (1 - x) * (b - a) / 2 + a;
			points[n-1-i] = (1 + x) * (b - a) / 2 + a;
		}
		if (n % 2 == 1)
			points[n/2] = (a+b)/2;
		return points;
	}

	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> sigmoid(SizeT n, Number a, Number b, Number steepness) {
		if (n == 1)
			return {(a+b)/2};
		Container<Number> points(n);
		for (SizeT i = 0; i < n; ++i) {
			const Number x = static_cast<Number>(i) / (static_cast<Number>(n) - 1);
			const Number sigmoidValue = 1.0 / (1.0 + exp(-steepness * (x - 0.5)));
			points[i] = a + (b - a) * sigmoidValue;
		}
		return points;
	}

	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> sigmoid_stretched(SizeT n, Number a, Number b, Number steepness) {
		if (n == 0)
			return {};
		if (n == 1)
			return {(a+b)/2};
		Container<Number> points(n);
		const Number stretch_factor = 1 - 2 / (1 + std::exp(0.5 * steepness));
		for (SizeT i = 1; i < n - 1; ++i) {
			const Number x = static_cast<double>(i) / (static_cast<Number>(n) - 1);
			const Number sigmoid = 1.0 / (1.0 + std::exp(-steepness * (x - 0.5)));
			const Number stretched = (sigmoid - 1.0 / (1.0 + std::exp(0.5 * steepness))) / stretch_factor;
			points[i] = a + (b - a) * stretched;
		}
		points[0] = a;
		points[n-1] = b;
		return points;
	}

	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> erf(SizeT n, Number a, Number b, Number steepness) {
		if (n == 0)
			return {};
		if (n == 1)
			return {(a+b)/2};
		Container<Number> points(n);
		for (SizeT i = 0; i < n; ++i) {
			const Number x = static_cast<Number>(i) / (static_cast<Number>(n) - 1);
			const Number erf_value = std::erf(steepness * (x - 0.5));
			points[i] = a + (b - a) * (1 + erf_value) / 2;
		}
		return points;
	}

	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> erf_stretched(SizeT n, Number a, Number b, Number steepness) {
		Container<Number> points = erf(n, a, b, steepness);
		stretch(points, a, b);
		return points;
	}

	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> of_type(Type type, SizeT n, Number a, Number b, Number parameter = 0) {
		switch (type) {
		case Type::CHEBYSHEV:
			return chebyshev(n, a, b);
		case Type::CHEBYSHEV_STRETCHED:
			return chebyshev_stretched(n, a, b);
		case Type::CHEBYSHEV_ELLIPSE:
			return chebyshev_ellipse(n, a, b, parameter);
		case Type::CHEBYSHEV_ELLIPSE_STRETCHED:
			return chebyshev_ellipse_stretched(n, a, b, parameter);
		case Type::CIRCLE_PROJECTION:
			return circle_proj(n, a, b);
		case Type::ELLIPSE_PROJECTION:
			return ellipse_proj(n, a, b, parameter);
		case Type::SIGMOID:
			return sigmoid(n, a, b, parameter);
		case Type::SIGMOID_STRETCHED:
			return sigmoid_stretched(n, a, b, parameter);
		case Type::ERF:
			return erf(n, a, b, parameter);
		case Type::ERF_STRETCHED:
			return erf_stretched(n, a, b, parameter);
		case Type::UNIFORM:
		case Type::GENERAL:
		default:
			return uniform(n, a, b);
		}
	}
}