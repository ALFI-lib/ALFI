#pragma once

#ifndef ALFI_DEFAULT_NUMBER
#define ALFI_DEFAULT_NUMBER double
#endif

#ifndef ALFI_DEFAULT_CONTAINER
#include <vector>
#define ALFI_DEFAULT_CONTAINER std::vector<Number>
#endif

namespace alfi {
	using DefaultNumber = ALFI_DEFAULT_NUMBER;

	template <typename Number = DefaultNumber>
	using DefaultContainer = ALFI_DEFAULT_CONTAINER;
}