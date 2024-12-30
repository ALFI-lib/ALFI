#pragma once

#include <iostream>

#include "../config.h"

namespace alfi::util::arrays {
	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> add(const Container<Number>& container1, const Container<Number>& container2) {
		if (container1.size() != container2.size()) {
			std::cerr << "Error in function " << __FUNCTION__
					  << ": Vectors container1 (of size " << container1.size()
					  << ") and container2 (of size " << container2.size()
					  << ") are not the same size. Returning an empty array..." << std::endl;
			return {};
		}
		const auto n = container1.size();
		Container<Number> result(n);
		for (std::remove_const_t<decltype(n)> i = 0; i < n; ++i) {
			result[i] = container1[i] + container2[i];
		}
		return result;
	}

	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> sub(const Container<Number>& container1, const Container<Number>& container2) {
		if (container1.size() != container2.size()) {
			std::cerr << "Error in function " << __FUNCTION__
					  << ": Vectors container1 (of size " << container1.size()
					  << ") and container2 (of size " << container2.size()
					  << ") are not the same size. Returning an empty array..." << std::endl;
			return {};
		}
		const auto n = container1.size();
		Container<Number> result(n);
		for (std::remove_const_t<decltype(n)> i = 0; i < n; ++i) {
			result[i] = container1[i] - container2[i];
		}
		return result;
	}

	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> diff(const auto& container) {
		if (container.empty()) {
			return {};
		}
		Container<Number> result(container.size() - 1);
		for (SizeT i = 0; i < result.size(); ++i) {
			result[i] = container[i+1] - container[i];
		}
		return result;
	}

	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> mean(const Container<Number>& container1, const Container<Number>& container2) {
		if (container1.size() != container2.size()) {
			std::cerr << "Error in function " << __FUNCTION__
					  << ": Vectors container1 (of size " << container1.size()
					  << ") and container2 (of size " << container2.size()
					  << ") are not the same size. Returning an empty array..." << std::endl;
			return {};
		}
		const auto n = container1.size();
		Container<Number> result(n);
		for (std::remove_const_t<decltype(n)> i = 0; i < n; ++i) {
			result[i] = (container1[i] + container2[i]) / 2;
		}
		return result;
	}
}