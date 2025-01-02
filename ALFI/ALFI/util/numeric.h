#pragma once

#include "../config.h"

namespace alfi::util::numeric {
	template <typename Number = DefaultNumber>
	bool are_equal(Number a, Number b, Number epsilon = std::numeric_limits<Number>::epsilon()) {
		return std::abs(a - b) < epsilon || std::abs(a - b) < std::max(std::abs(a), std::abs(b)) * epsilon;
	}
}