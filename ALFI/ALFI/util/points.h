#pragma once

#include "../config.h"

namespace alfi::points {
	template <typename Number = DefaultNumber>
	void lin_map(auto& points, Number a, Number b, Number c, Number d) {
		const auto mid1 = (a + b) / 2;
		const auto mid2 = (c + d) / 2;
		const auto scale = (d - c) / (b - a);
		for (auto& point : points) {
			point = mid2 + scale * (point - mid1);
		}
	}

	template <typename Number = DefaultNumber>
	auto lin_mapped(const auto& points, Number a, Number b, Number c, Number d) {
		auto mapped_points = points;
		lin_map(mapped_points, a, b, c, d);
		return mapped_points;
	}

	template <typename Number = DefaultNumber>
	void stretch(auto& points, Number a, Number b) {
		if (points.empty()) {
			return;
		}
		if (points.size() == 1 || points.front() == points.back()) {
			std::fill(points.begin(), points.end(), (a + b) / 2);
			return;
		}
		lin_map(points, points.front(), points.back(), a, b);
		points.front() = a;
		points.back() = b;
	}

	template <typename Number = DefaultNumber>
	auto stretched(const auto& points, Number a, Number b) {
		auto stretched_points = points;
		stretch(stretched_points, a, b);
		return stretched_points;
	}
}