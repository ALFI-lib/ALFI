#pragma once

#include <functional>

namespace alfi::util {
	class SimpleScopeGuard {
	public:
		// ReSharper disable once CppNonExplicitConvertingConstructor
		SimpleScopeGuard(std::function<void()> on_exit) : _on_exit(std::move(on_exit)) {} // NOLINT(*-explicit-constructor)
		~SimpleScopeGuard() noexcept { _on_exit(); }
		SimpleScopeGuard(const SimpleScopeGuard&) = delete;
		SimpleScopeGuard& operator=(const SimpleScopeGuard&) = delete;
	private:
		std::function<void()> _on_exit;
	};

	template <typename Iterator, typename T>
	Iterator first_leq_or_begin(Iterator begin, Iterator end, const T& value) {
		auto iter = std::lower_bound(begin, end - 1, value);
		if (iter != begin && *iter > value) {
			--iter;
		}
		return iter;
	}
}