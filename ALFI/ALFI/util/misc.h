#pragma once

#include <algorithm>
#include <utility>

namespace alfi::util::misc {
	template <typename Function>
	class SimpleScopeGuard {
	public:
		// ReSharper disable once CppNonExplicitConvertingConstructor
		SimpleScopeGuard(Function on_exit) : _on_exit(std::move(on_exit)) {} // NOLINT(*-explicit-constructor)
		~SimpleScopeGuard() noexcept { _on_exit(); }
		SimpleScopeGuard(const SimpleScopeGuard&) = delete;
		SimpleScopeGuard& operator=(const SimpleScopeGuard&) = delete;
	private:
		Function _on_exit;
	};

	template <typename... Ts>
	struct overload : Ts... {
		using Ts::operator()...;
	};
	template <typename... Ts>
	overload(Ts...) -> overload<Ts...>;

	template <typename Iterator, typename T>
	Iterator first_leq_or_begin(Iterator begin, Iterator end, const T& value) {
		auto iter = std::lower_bound(begin, end - 1, value);
		if (iter != begin && *iter > value) {
			--iter;
		}
		return iter;
	}
}