#pragma once

#ifndef ALFI_DEFAULT_NUMBER
#define ALFI_DEFAULT_NUMBER double
#endif

#ifndef ALFI_DEFAULT_CONTAINER
#include <vector>
#define ALFI_DEFAULT_CONTAINER std::vector
#endif

#ifndef ALFI_SIZE_TYPE
#define ALFI_SIZE_TYPE size_t
#endif

namespace alfi {
	using DefaultNumber = ALFI_DEFAULT_NUMBER;

	template <typename Number = DefaultNumber, typename... Ts>
	using DefaultContainer = ALFI_DEFAULT_CONTAINER<Number, Ts...>;

	using SizeT = ALFI_SIZE_TYPE;

	namespace traits {
		template <typename T, typename = void>
		struct has_size : std::false_type {};
		template <typename T>
		struct has_size<T, std::void_t<decltype(std::declval<T>().size())>> : std::true_type {};
	}
}